#include "bhc/bhc_c_api.h"
#include "bhc/bhc.hpp"
#include "common.hpp"  // 为了 GetFieldAddr

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

// 使用 Bellhop C++ 命名空间
using namespace bhc;

// 线程局部存储，用于保存最后一次错误信息
static thread_local std::string g_last_error_message;

// 内部使用的辅助函数
namespace {

// 将 C 结构体 BhcSSP 转换为 Bellhop 内部参数
void setup_ssp(bhcParams<false>& params, const BhcSSP* ssp)
{
    if (ssp == nullptr || ssp->n_points < 2) {
        throw std::invalid_argument("SSP 参数非法：至少需要 2 个点");
    }
    if (ssp->z == nullptr || ssp->c == nullptr || ssp->rho == nullptr) {
        throw std::invalid_argument("SSP 参数非法：z/c/rho 指针不能为空");
    }

    // Bellhop 内部使用 1D SSP 数组存储，并在预处理阶段根据类型生成插值系数
    params.ssp->Type = static_cast<char>(ssp->type);
    params.ssp->NPts = ssp->n_points;
    params.ssp->dirty = true; // 标记为需要预处理

    for (int i = 0; i < ssp->n_points; ++i) {
        params.ssp->z[i] = static_cast<real>(ssp->z[i]);
        params.ssp->alphaR[i] = static_cast<real>(ssp->c[i]); // 声速存储在 alphaR
        params.ssp->rho[i] = static_cast<real>(ssp->rho[i]);
        params.ssp->alphaI[i] = 0; // 默认为0
        if (ssp->alpha != nullptr) {
            // 衰减系数的单位转换逻辑可以根据需要在此处添加
            params.ssp->alphaI[i] = static_cast<real>(ssp->alpha[i]);
        }
    }
}

// 将 C 结构体 BhcBathymetry 转换为 Bellhop 内部参数
void setup_bathymetry(bhcParams<false>& params, const BhcBathymetry* bty)
{
    if (bty == nullptr || bty->n_points <= 0) {
        // 如果没有提供地形，则使用平坦海底（深度为SSP的最深点）
        params.Bdry->Bot.hs.Depth = params.ssp->z[params.ssp->NPts - 1];
        return;
    }
    if (bty->z == nullptr) {
        throw std::invalid_argument("BTY 参数非法：z 指针不能为空");
    }

    // 当前实现简化为取地形中的最大深度作为平坦海底深度
    double max_bty_depth = 0.0;
    for (int i = 0; i < bty->n_points; ++i) {
        max_bty_depth = std::max(max_bty_depth, bty->z[i]);
    }

    double max_ssp_depth = params.ssp->z[params.ssp->NPts - 1];
    params.Bdry->Bot.hs.Depth = static_cast<real>(std::min(max_ssp_depth, max_bty_depth));
}

// 将 C 结构体 BhcConfig 转换为 Bellhop 内部参数
void setup_config(bhcParams<false>& params, const BhcConfig* config)
{
    if (config == nullptr) {
        throw std::invalid_argument("Config 不能为空");
    }

    // 标题
    std::memset(params.Title, 0, sizeof(params.Title));
    if (config->title != nullptr) {
        strncpy(params.Title, config->title, sizeof(params.Title) - 1);
    } else {
        strncpy(params.Title, "BHC C API - TL Calculation", sizeof(params.Title) - 1);
    }

    // RunType 字符串（7个字符）
    std::memset(params.Beam->RunType, ' ', 7);
    params.Beam->RunType[0] = static_cast<char>(config->run_type);
    params.Beam->RunType[1] = static_cast<char>(config->beam.beam_type);
    params.Beam->RunType[3] = 'R'; // 点源
    params.Beam->RunType[4] = 'R'; // 规则网格
    params.Beam->RunType[5] = ' '; // 2D

    // 频率
    params.freqinfo->freq0 = static_cast<real>(config->frequency);

    // 源深度
    extsetup_sz<false>(params, 1);
    params.Pos->Sz[0] = static_cast<float>(config->src_depth);

    // 接收器网格
    if (config->n_ranges < 1 || config->n_depths < 1) {
        throw std::invalid_argument("接收器网格非法：n_ranges/n_depths 必须 >= 1");
    }
    extsetup_rcvrranges<false>(params, config->n_ranges);
    extsetup_rcvrdepths<false>(params, config->n_depths);

    if (config->n_ranges == 1) {
        params.Pos->Rr[0] = static_cast<float>(config->range_min);
    } else {
        double dr = (config->range_max - config->range_min) / double(config->n_ranges - 1);
        for (int i = 0; i < config->n_ranges; ++i) {
            params.Pos->Rr[i] = static_cast<float>(config->range_min + double(i) * dr);
        }
    }

    if (config->n_depths == 1) {
        params.Pos->Rz[0] = static_cast<float>(config->depth_min);
    } else {
        double dz = (config->depth_max - config->depth_min) / double(config->n_depths - 1);
        for (int i = 0; i < config->n_depths; ++i) {
            params.Pos->Rz[i] = static_cast<float>(config->depth_min + double(i) * dz);
        }
    }

    // 发射角
    if (config->beam.n_beams < 1) {
        throw std::invalid_argument("beam.n_beams 必须 >= 1");
    }
    extsetup_rayelevations<false>(params, config->beam.n_beams);
    params.Angles->alpha.inDegrees = true;
    if (config->beam.n_beams == 1) {
        params.Angles->alpha.angles[0] = static_cast<real>(config->beam.angle_min);
    } else {
        double da = (config->beam.angle_max - config->beam.angle_min) / double(config->beam.n_beams - 1);
        for (int i = 0; i < config->beam.n_beams; ++i) {
            params.Angles->alpha.angles[i] = static_cast<real>(config->beam.angle_min + double(i) * da);
        }
    }

    // 边界条件
    params.Bdry->Top.hs.Depth = 0.0;
    params.Bdry->Top.hs.bc = static_cast<char>(config->boundary.surface_type);
    params.Bdry->Bot.hs.bc = static_cast<char>(config->boundary.bottom_type);

    if (config->boundary.bottom_type == BHC_BOUNDARY_ACOUSTIC) {
        params.Bdry->Bot.hs.rho = static_cast<real>(config->boundary.bottom_rho);
        params.Bdry->Bot.hs.alphaR = static_cast<real>(config->boundary.bottom_cp);
        params.Bdry->Bot.hs.betaR = static_cast<real>(config->boundary.bottom_cs);
        params.Bdry->Bot.hs.alphaI = static_cast<real>(config->boundary.bottom_ap);
        params.Bdry->Bot.hs.betaI = static_cast<real>(config->boundary.bottom_as);
    }
}

// 从 Bellhop 输出中提取 TL 结果
void extract_tl_results(const bhcOutputs<false, false>& outputs, const bhcParams<false>& params, BhcTLResult* result)
{
    int n_ranges = params.Pos->NRr;
    int n_depths = params.Pos->NRz;

    result->n_ranges = n_ranges;
    result->n_depths = n_depths;
    result->n_freqs = 1;
    result->n_bearings = 1;
    result->tl = new float[n_ranges * n_depths];

    // 2D情况下，源索引和方位角索引都为0
    for (int ir = 0; ir < n_ranges; ++ir) {
        for (int iz = 0; iz < n_depths; ++iz) {
            // 计算一维索引
            size_t idx = GetFieldAddr(0, 0, 0, 0, iz, ir, params.Pos);
            
            // 获取复声压并计算 TL (dB)
            cpxf p = outputs.uAllSources[idx];
            float amp = std::sqrt(p.real() * p.real() + p.imag() * p.imag());
            result->tl[iz * n_ranges + ir] = (amp > 0) ? -20.0f * log10f(amp) : 200.0f;
        }
    }
}

} // 匿名命名空间结束

// --- C API 公开函数实现 ---

BHC_C_API BhcTLResult* bhc_compute_tl(const BhcConfig* config, const BhcSSP* ssp, const BhcBathymetry* bty)
{
    g_last_error_message.clear();
    try {
        // 1. 初始化 Bellhop
        bhcInit init{}; // 使用 C++11 的值初始化，确保所有成员都被正确初始化（包括默认值）
        init.FileRoot = nullptr;
        init.prtCallback = [](const char* msg) {
            if (msg) g_last_error_message.append(msg);
        };
        init.outputCallback = [](const char* msg) { 
            if (msg) g_last_error_message.append(msg); 
        };

        bhcParams<false> params;
        bhcOutputs<false, false> outputs;

        // 初始化（不读文件，但会填充各模块默认值）
        bool setup_ok = bhc::setup_nofile<false, false>(init, params, outputs);
        if(!setup_ok) {
            if(g_last_error_message.empty()) {
                g_last_error_message = "初始化失败：bhc::setup_nofile 返回 false";
            }
            return nullptr;
        }

        // 2. 将 C 结构体转换为 Bellhop 内部参数
        setup_config(params, config);
        setup_ssp(params, ssp);
        setup_bathymetry(params, bty);

        // 3. 验证参数并回显（可选，但有助于调试）
        bool echo_ok = bhc::echo<false>(params);
        if (!echo_ok) {
            bhc::finalize<false, false>(params, outputs);
            if (g_last_error_message.empty()) {
                g_last_error_message = "参数验证失败（未收到更详细的错误信息）";
            }
            return nullptr;
        }

        // 4. 运行计算
        bool run_ok = bhc::run<false, false>(params, outputs);
        if (!run_ok) {
            bhc::finalize<false, false>(params, outputs);
            g_last_error_message = "Bellhop run 失败";
            return nullptr;
        }

        // 5. 提取结果
        BhcTLResult* result = new BhcTLResult();
        std::memset(result, 0, sizeof(BhcTLResult));
        extract_tl_results(outputs, params, result);

        // 6. 清理 Bellhop 资源
        bhc::finalize<false, false>(params, outputs);

        return result;

    } catch (const std::exception& e) {
        g_last_error_message = e.what();
        return nullptr;
    } catch (...) {
        g_last_error_message = "发生了未知异常";
        return nullptr;
    }
}

BHC_C_API void bhc_free_tl_result(BhcTLResult* result)
{
    if (result) {
        delete[] result->tl;
        delete result;
    }
}

BHC_C_API const char* bhc_get_last_error(void)
{
    return g_last_error_message.empty() ? nullptr : g_last_error_message.c_str();
}

BHC_C_API const char* bhc_get_version(void)
{
    return "Bellhop C API v1.0";
}

BHC_C_API int bhc_initialize(void) { return 0; }
BHC_C_API void bhc_cleanup(void) {}
