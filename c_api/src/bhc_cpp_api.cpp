#include "bhc/bhc_cpp_api.hpp"
#include "bhc/bhc.hpp"
#include "common.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>

namespace bhc_cpp {

// ... (辅助函数: pad_right, set_title, set_runtype, set_beam, set_freq0, set_positions,
// set_angles) ...

static inline std::string pad_right(std::string s, size_t n)
{
    if(s.size() >= n) {
        s.resize(n);
        return s;
    }
    s.append(n - s.size(), ' ');
    return s;
}
static inline void set_title(bhc::bhcParams<false> &params, const Title &t)
{
    std::memset(params.Title, 0, sizeof(params.Title));
    std::string s = t.text;
    if(s.empty()) s = "bhc_cpp";
    if(s.size() >= sizeof(params.Title)) s.resize(sizeof(params.Title) - 1);
    std::memcpy(params.Title, s.data(), s.size());
}
static inline void set_runtype(bhc::bhcParams<false> &params, const RunType &rt)
{
    char out[7];
    for(int i = 0; i < 7; ++i) out[i] = ' ';
    if(rt.use_raw) {
        std::string s = pad_right(rt.raw, 7);
        std::memcpy(out, s.data(), 7);
    } else {
        out[0] = static_cast<char>(rt.e.tl);
        out[1] = static_cast<char>(rt.e.infl);
        out[3] = static_cast<char>(rt.e.source);
        out[4] = static_cast<char>(rt.e.grid);
        out[5] = static_cast<char>(rt.e.dim);
    }
    std::memcpy(params.Beam->RunType, out, 7);
}
static inline void set_beam(bhc::bhcParams<false> &params, const Beam &b)
{
    std::string t = pad_right(b.type, 4);
    std::memcpy(params.Beam->Type, t.data(), 4);
    params.Beam->Box.x     = static_cast<bhc::real>(b.box_x_m);
    params.Beam->Box.y     = static_cast<bhc::real>(b.box_z_m);
    params.Beam->rangeInKm = false;
    if(b.deltas_m > 0) {
        params.Beam->autoDeltas = false;
        params.Beam->deltas     = static_cast<bhc::real>(b.deltas_m);
    } else {
        params.Beam->autoDeltas = true;
        params.Beam->deltas     = 0;
    }
    params.Beam->epsMultiplier = static_cast<bhc::real>(b.eps_multiplier);
}
static inline void set_freq0(bhc::bhcParams<false> &params, const Freq0 &f)
{
    params.freqinfo->freq0 = static_cast<bhc::real>(f.hz);
}
static inline void set_positions(bhc::bhcParams<false> &params, const Positions2D &p)
{
    if(p.sz_m.empty() || p.rr_m.empty() || p.rz_m.empty())
        throw Error("Positions2D 向量不能为空");
    bhc::extsetup_sz<false>(params, p.sz_m.size());
    for(size_t i = 0; i < p.sz_m.size(); ++i) params.Pos->Sz[i] = p.sz_m[i];
    bhc::extsetup_rcvrranges<false>(params, p.rr_m.size());
    for(size_t i = 0; i < p.rr_m.size(); ++i) params.Pos->Rr[i] = p.rr_m[i];
    bhc::extsetup_rcvrdepths<false>(params, p.rz_m.size());
    for(size_t i = 0; i < p.rz_m.size(); ++i) params.Pos->Rz[i] = p.rz_m[i];
    params.Pos->RrInKm = p.rr_in_km;
}
static inline void set_angles(
    bhc::bhcParams<false> &params, const Angles &a, int n_rays, double start_angle_deg,
    double end_angle_deg)
{
    // 兼容两种输入方式：
    // 1) 显式给 angles.alpha（优先）
    // 2) angles.alpha 为空时，用 (n_rays, start_angle_deg, end_angle_deg) 生成
    std::vector<double> alpha;
    bool in_degrees = true;

    if(!a.alpha.empty()) {
        alpha = a.alpha;
        in_degrees = a.alpha_in_degrees;
    } else {
        if(n_rays <= 0) throw Error("Angles.alpha 为空且 n_rays<=0，无法生成射线角");
        alpha.resize(static_cast<size_t>(n_rays));
        for(int i = 0; i < n_rays; ++i) {
            double t = (n_rays == 1) ? 0.0 : double(i) / double(n_rays - 1);
            alpha[static_cast<size_t>(i)] = start_angle_deg + (end_angle_deg - start_angle_deg) * t;
        }
        in_degrees = true; // env 的 START/END ANGLES 默认是度
    }

    bhc::extsetup_rayelevations<false>(params, static_cast<int32_t>(alpha.size()));
    params.Angles->alpha.inDegrees = in_degrees;
    for(size_t i = 0; i < alpha.size(); ++i) {
        params.Angles->alpha.angles[i] = alpha[i];
    }
}

// 设置 SSP 深度网格和 1D 属性
static inline void set_ssp_depth_grid(bhc::bhcParams<false> &params, const SSP1D &ssp)
{
    if(ssp.pts.size() < 2) throw Error("SSP1D.pts 至少需要 2 个点");

    params.ssp->NPts = static_cast<int32_t>(ssp.pts.size());
    params.ssp->Nz   = params.ssp->NPts;

    std::string au           = pad_right(ssp.atten_unit, 2);
    params.ssp->AttenUnit[0] = au[0];
    params.ssp->AttenUnit[1] = au[1];

    for(int32_t i = 0; i < params.ssp->NPts; ++i) {
        const auto &p         = ssp.pts[i];
        params.ssp->z[i]      = p.z_m;
        params.ssp->alphaR[i] = p.alphaR_mps;
        params.ssp->betaR[i]  = p.betaR_mps;
        params.ssp->rho[i]    = p.rho_gcm3;
        params.ssp->alphaI[i] = p.alphaI;
        params.ssp->betaI[i]  = p.betaI;
    }

    params.ssp->dirty         = true;
    params.Bdry->Top.hs.Depth = params.ssp->z[0];
    params.Bdry->Bot.hs.Depth = params.ssp->z[params.ssp->NPts - 1];
}

// 设置 2D Range-Dependent SSP (Type 'Q')
static inline void set_ssp_quad(
    bhc::bhcParams<false> &params, const SSP1D &ssp_base, const SSPQuad &ssp_quad)
{
    const int32_t n_depths = ssp_base.pts.size();
    const int32_t n_ranges = ssp_quad.ranges.size();

    if(n_depths < 2 || n_ranges < 2) {
        throw Error("SSP Type 'Q' 至少需要 2个深度和2个距离点");
    }
    if(ssp_quad.ssp_matrix.size() != static_cast<size_t>(n_depths) * n_ranges) {
        throw Error("SSPQuad.ssp_matrix 尺寸与深度/距离点数不匹配");
    }

    // 调用底层 API 分配内存
    bhc::extsetup_ssp_quad(params, n_depths, n_ranges);

    // 填充深度和距离向量
    for(int i = 0; i < n_depths; ++i) params.ssp->z[i] = ssp_base.pts[i].z_m;
    for(int i = 0; i < n_ranges; ++i) params.ssp->Seg.r[i] = ssp_quad.ranges[i];

    // 填充声速矩阵 (cMat 是平铺的，布局与我们的 ssp_matrix 一致)
    for(size_t i = 0; i < ssp_quad.ssp_matrix.size(); ++i) {
        params.ssp->cMat[i] = ssp_quad.ssp_matrix[i];
    }

    params.ssp->rangeInKm = ssp_quad.ranges_in_km;
    params.ssp->dirty     = true;
}

// ... (set_boundaries_2d, write_boundary_curve_2d) ...
static inline void write_boundary_curve_2d(
    bhc::bhcParams<false> &params, bhc::BdryInfoTopBot<false> &dst, const Boundary2D &src,
    bool is_top)
{
    if(src.r.size() != src.z.size() || src.r.size() < 2)
        throw Error("Boundary2D 点列非法");
    const char t0 = static_cast<char>(src.interp);
    if(t0 != 'C' && t0 != 'L') throw Error("2D Boundary.interp 必须是 'C' 或 'L'");
    dst.type[0]        = t0;
    dst.type[1]        = static_cast<char>(src.flag);
    dst.rangeInKm      = src.range_in_km;
    dst.dirty          = true;
    const int32_t N    = src.r.size();
    const bool ext     = src.extend_to_infinity;
    const int32_t NPts = ext ? (N + 2) : N;
    if(is_top)
        bhc::extsetup_altimetry<false>(params, NPts);
    else
        bhc::extsetup_bathymetry<false>(params, NPts);
    auto write_point = [&](int32_t i, double rr, double zz) {
        dst.bd[i].x = bhc::vec2(rr, zz);
    };
    int32_t offset = ext ? 1 : 0;
    if(ext) write_point(0, -src.extend_left, src.z.front());
    for(int32_t i = 0; i < N; ++i) write_point(i + offset, src.r[i], src.z[i]);
    if(ext) write_point(NPts - 1, src.extend_right, src.z.back());
}
static inline void set_halfspace(bhc::BdryPtSmall &dst, const Halfspace &src)
{
    dst.hs.bc     = static_cast<char>(src.bc);
    dst.hs.alphaR = static_cast<bhc::real>(src.alphaR_mps);
    dst.hs.betaR  = static_cast<bhc::real>(src.betaR_mps);
    dst.hs.rho    = static_cast<bhc::real>(src.rho_gcm3);
    dst.hs.alphaI = static_cast<bhc::real>(src.alphaI);
    dst.hs.betaI  = static_cast<bhc::real>(src.betaI);

    // Opt[6]
    std::string opt = pad_right(src.opt, 6);
    std::memcpy(dst.hs.Opt, opt.data(), 6);

    // Sigma（echo but usually not used）
    dst.hsx.Sigma = static_cast<bhc::real>(src.sigma);
}

static inline void set_boundaries_2d(bhc::bhcParams<false> &params, const Boundaries2D &b)
{
    set_halfspace(params.Bdry->Top, b.top);
    set_halfspace(params.Bdry->Bot, b.bot);

    if(b.top_curve.r.empty() || b.bot_curve.r.empty())
        throw Error("Boundary curves 必须提供");
    write_boundary_curve_2d(params, params.bdinfo->top, b.top_curve, true);
    write_boundary_curve_2d(params, params.bdinfo->bot, b.bot_curve, false);
}

BHC_CPP_API TLResult2D compute_tl_2d(const Input2D &in)
{
    bhc::bhcInit init{};
    init.FileRoot    = nullptr;
    init.prtCallback = [](const char *m) {
        if(m) std::cerr << m;
    };
    init.outputCallback = [](const char *m) {
        if(m) std::cerr << m;
    };

    bhc::bhcParams<false> params;
    bhc::bhcOutputs<false, false> outputs;

    if(!bhc::setup_nofile<false, false>(init, params, outputs)) {
        throw Error("bhc::setup_nofile 失败");
    }

        // --- 强制 no-file 模式，杜绝任何文件读取 ---
    // 1) 将 Top/Bottom 边界 Opt[6] 全清空，然后把 TopOpt 第 0 位设为 'V'
    //    这样 TopOpt='V '，底层会直接走 Vacuum 分支，不会尝试读 SSP 文件。
    std::memset(params.Bdry->Top.hs.Opt, ' ', sizeof(params.Bdry->Top.hs.Opt));
    std::memset(params.Bdry->Bot.hs.Opt, ' ', sizeof(params.Bdry->Bot.hs.Opt));
    params.Bdry->Top.hs.Opt[0] = 'V';

    // 2) 标记 altimetry/bathymetry 已修改，确保 preprocess 重新计算但不读文件。
    params.bdinfo->top.dirty = true;
    params.bdinfo->bot.dirty = true;

    // --- 参数转换 ---
    set_title(params, in.title);
    set_freq0(params, in.freq0);
    // 允许外部设置 top boundary 的 Opt，但为了 nofile 模式安全，不允许 Opt[0] 触发从文件读取 SSP。
    // 如果用户需要 Q/H SSP，必须通过内存接口提供（ssp_quad/ssp_matrix）。
    // 注：TopOpt::Validate 已在底层对 nofile 模式跳过 .ssp 文件存在性检查。
    set_runtype(params, in.run_type);
    set_beam(params, in.beam);
    set_positions(params, in.pos);
    set_angles(params, in.angles, in.n_rays, in.start_angle_deg, in.end_angle_deg);

    // 根据 options 决定 SSP 类型
    bool is_quad_ssp = (in.options.find('Q') != std::string::npos);
    if(is_quad_ssp) {
        params.ssp->Type = 'Q';
        set_ssp_depth_grid(params, in.ssp);        // 先用 1D SSP 设置深度网格
        set_ssp_quad(params, in.ssp, in.ssp_quad); // 再用 2D SSP 填充数据
    } else {
        params.ssp->Type = static_cast<char>(in.ssp.type);
        set_ssp_depth_grid(params, in.ssp);
    }

    set_boundaries_2d(params, in.boundaries);

    // --- 运行与结果提取 ---
    if(!bhc::echo<false>(params)) {
        bhc::finalize<false, false>(params, outputs);
        throw Error("bhc::echo 失败（请查看控制台输出）");
    }

    if(!bhc::run<false, false>(params, outputs)) {
        bhc::finalize<false, false>(params, outputs);
        throw Error("bhc::run 失败（请查看控制台输出）");
    }

    TLResult2D out;
    out.width  = params.Pos->NRr;
    out.height = params.Pos->NRz;
    out.tl_db.resize(out.width * out.height);

    for(int ir = 0; ir < out.width; ++ir) {
        for(int iz = 0; iz < out.height; ++iz) {
            size_t idx = bhc::GetFieldAddr(0, 0, 0, 0, iz, ir, params.Pos);
            auto p     = outputs.uAllSources[idx];
            float amp  = std::sqrt(p.real() * p.real() + p.imag() * p.imag());
            out.tl_db[iz * out.width + ir] = (amp > 0) ? -20.0f * std::log10(amp)
                                                       : 200.0f;
        }
    }

    bhc::finalize<false, false>(params, outputs);
    return out;
}
} // namespace bhc_cpp
