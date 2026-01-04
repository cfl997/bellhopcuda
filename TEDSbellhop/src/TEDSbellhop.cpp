#include "TEDSbellhop/TEDSbellhop.hpp"

#include "bhc/bhc.hpp"   // 内部使用 bellhopcxx/bellhopcuda C++ API
#include "common.hpp"    // 与 c_api 保持一致（工程内已有）

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

namespace tedsbellhop {

// 进度：按线程保存“当前正在运行的 params 指针”，以及一个缓存的百分比。
// 约定：当没有正在运行的计算时，g_progress_params == nullptr，get_percent_progress() 返回 -1。
static thread_local bhc::bhcParams<false>* g_progress_params = nullptr;
static thread_local int g_progress_cache = -1;

int get_percent_progress()
{
    if(!g_progress_params) return -1;
    // bhc::get_percent_progress 是线程安全的，返回 0~100
    g_progress_cache = bhc::get_percent_progress<false>(*g_progress_params);
    return g_progress_cache;
}

static inline std::string pad_right(std::string s, size_t n)
{
    if(s.size() >= n) {
        s.resize(n);
        return s;
    }
    s.append(n - s.size(), ' ');
    return s;
}

static inline void set_title(bhc::bhcParams<false> &params, const std::string &title)
{
    std::memset(params.Title, 0, sizeof(params.Title));
    std::string s = title;
    if(s.empty()) s = "TEDSbellhop";
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
        out[0] = static_cast<char>(rt.tl);
        out[1] = static_cast<char>(rt.infl);
        out[3] = static_cast<char>(rt.source);
        out[4] = static_cast<char>(rt.grid);
        out[5] = static_cast<char>(rt.dim);
    }

    std::memcpy(params.Beam->RunType, out, 7);
}

static inline void set_beam(bhc::bhcParams<false> &params, const Beam &b)
{
    std::string t = pad_right(b.type, 4);
    std::memcpy(params.Beam->Type, t.data(), 4);

    params.Beam->Box.x = static_cast<bhc::real>(b.box_x_m);
    params.Beam->Box.y = static_cast<bhc::real>(b.box_z_m);

    // 我们对外全部用 m，内部用 rangeInKm 标记控制范围单位
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

static inline void set_freq0(bhc::bhcParams<false> &params, double hz)
{
    params.freqinfo->freq0 = static_cast<bhc::real>(hz);
}

static inline void set_positions(bhc::bhcParams<false> &params, const Positions2D &p)
{
    if(p.source_depths_m.empty() || p.receiver_ranges.empty() || p.receiver_depths_m.empty())
        throw Error("Positions2D 向量不能为空 (source_depths_m/receiver_ranges/receiver_depths_m)");

    bhc::extsetup_sz<false>(params, p.source_depths_m.size());
    for(size_t i = 0; i < p.source_depths_m.size(); ++i) params.Pos->Sz[i] = p.source_depths_m[i];

    bhc::extsetup_rcvrranges<false>(params, p.receiver_ranges.size());
    for(size_t i = 0; i < p.receiver_ranges.size(); ++i) params.Pos->Rr[i] = p.receiver_ranges[i];

    bhc::extsetup_rcvrdepths<false>(params, p.receiver_depths_m.size());
    for(size_t i = 0; i < p.receiver_depths_m.size(); ++i) params.Pos->Rz[i] = p.receiver_depths_m[i];

    params.Pos->RrInKm = p.receiver_ranges_in_km;
}

static inline void set_angles(
    bhc::bhcParams<false> &params, const Angles &a, int n_rays, double start_angle_deg,
    double end_angle_deg)
{
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
        in_degrees = true;
    }

    bhc::extsetup_rayelevations<false>(params, static_cast<int32_t>(alpha.size()));
    params.Angles->alpha.inDegrees = in_degrees;
    for(size_t i = 0; i < alpha.size(); ++i) {
        params.Angles->alpha.angles[i] = alpha[i];
    }
}

static inline void set_ssp_depth_grid(bhc::bhcParams<false> &params, const SSP1D &ssp)
{
    if(ssp.points.size() < 2) throw Error("SSP1D.points 至少需要 2 个点");

    params.ssp->NPts = static_cast<int32_t>(ssp.points.size());
    params.ssp->Nz   = params.ssp->NPts;

    std::string au           = pad_right(ssp.atten_unit, 2);
    params.ssp->AttenUnit[0] = au[0];
    params.ssp->AttenUnit[1] = au[1];

    for(int32_t i = 0; i < params.ssp->NPts; ++i) {
        const auto &p = ssp.points[static_cast<size_t>(i)];
        params.ssp->z[i]      = p.z_m;
        params.ssp->alphaR[i] = p.c_mps;
        params.ssp->betaR[i]  = p.cs_mps;
        params.ssp->rho[i]    = p.rho_gcm3;
        params.ssp->alphaI[i] = p.atten_c;
        params.ssp->betaI[i]  = p.atten_s;
    }

    params.ssp->dirty         = true;
    params.Bdry->Top.hs.Depth = params.ssp->z[0];
    params.Bdry->Bot.hs.Depth = params.ssp->z[params.ssp->NPts - 1];
}

static inline void set_ssp_quad(
    bhc::bhcParams<false> &params, const SSP1D &ssp_base, const SSPQuad &ssp_quad)
{
    const int32_t n_depths = static_cast<int32_t>(ssp_base.points.size());
    const int32_t n_ranges = static_cast<int32_t>(ssp_quad.ranges.size());

    if(n_depths < 2 || n_ranges < 2) throw Error("SSPQuad 至少需要 2 个深度点和 2 个距离点");

    if(ssp_quad.c_mat.size() != static_cast<size_t>(n_depths) * static_cast<size_t>(n_ranges))
        throw Error("SSPQuad.c_mat 尺寸与深度/距离点数不匹配");

    bhc::extsetup_ssp_quad(params, n_depths, n_ranges);

    for(int i = 0; i < n_depths; ++i) params.ssp->z[i] = ssp_base.points[static_cast<size_t>(i)].z_m;
    for(int i = 0; i < n_ranges; ++i) params.ssp->Seg.r[i] = ssp_quad.ranges[static_cast<size_t>(i)];

    // 布局与 c_api 当前实现一致：直接平铺复制
    for(size_t i = 0; i < ssp_quad.c_mat.size(); ++i) params.ssp->cMat[i] = ssp_quad.c_mat[i];

    params.ssp->rangeInKm = ssp_quad.ranges_in_km;
    params.ssp->dirty     = true;
}

static inline void write_boundary_curve_2d(
    bhc::bhcParams<false> &params, bhc::BdryInfoTopBot<false> &dst, const BoundaryCurve2D &src,
    bool is_top)
{
    if(src.r.size() != src.z.size() || src.r.size() < 2)
        throw Error("BoundaryCurve2D 点列非法：r/z 长度不一致或点数<2");

    const char t0 = static_cast<char>(src.interp);
    if(t0 != 'C' && t0 != 'L') throw Error("BoundaryCurve2D.interp 必须是 'C' 或 'L'");

    dst.type[0]   = t0;
    dst.type[1]   = static_cast<char>(src.flag);
    dst.rangeInKm = src.range_in_km;
    dst.dirty     = true;

    const int32_t N = static_cast<int32_t>(src.r.size());
    const bool ext  = src.extend_to_infinity;
    const int32_t NPts = ext ? (N + 2) : N;

    if(is_top)
        bhc::extsetup_altimetry<false>(params, NPts);
    else
        bhc::extsetup_bathymetry<false>(params, NPts);

    auto write_point = [&](int32_t i, double rr, double zz) {
        dst.bd[i].x = bhc::vec2(rr, zz);
    };

    const int32_t offset = ext ? 1 : 0;
    if(ext) write_point(0, -src.extend_left, src.z.front());
    for(int32_t i = 0; i < N; ++i) write_point(i + offset, src.r[static_cast<size_t>(i)], src.z[static_cast<size_t>(i)]);
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

    std::string opt = pad_right(src.opt, 6);
    std::memcpy(dst.hs.Opt, opt.data(), 6);

    dst.hsx.Sigma = static_cast<bhc::real>(src.sigma);
}

static inline void set_boundaries_2d(bhc::bhcParams<false> &params, const Boundaries2D &b)
{
    set_halfspace(params.Bdry->Top, b.top);
    set_halfspace(params.Bdry->Bot, b.bottom);

    if(b.top_curve.r.empty() || b.bottom_curve.r.empty())
        throw Error("Boundaries2D 必须提供 top_curve/bottom_curve");

    write_boundary_curve_2d(params, params.bdinfo->top, b.top_curve, true);
    write_boundary_curve_2d(params, params.bdinfo->bot, b.bottom_curve, false);
}

void validate(const TL2DInput& in)
{
    if(in.pos.source_depths_m.empty()) throw Error("pos.source_depths_m 不能为空");
    if(in.pos.receiver_ranges.empty()) throw Error("pos.receiver_ranges 不能为空");
    if(in.pos.receiver_depths_m.empty()) throw Error("pos.receiver_depths_m 不能为空");

    if(in.ssp.points.size() < 2) throw Error("ssp.points 至少需要2个点");
    for(size_t i = 1; i < in.ssp.points.size(); ++i) {
        if(!(in.ssp.points[i].z_m > in.ssp.points[i-1].z_m))
            throw Error("ssp.points 的 z_m 必须严格单调递增");
    }

    const bool is_quad = (in.options.find('Q') != std::string::npos);
    if(is_quad) {
        if(in.ssp.type != SSPType::Quad)
            throw Error("options 包含 'Q' 时，ssp.type 必须为 SSPType::Quad");
        if(in.ssp_quad.ranges.size() < 2)
            throw Error("ssp_quad.ranges 至少需要2个点");
        const size_t n_depth = in.ssp.points.size();
        const size_t n_range = in.ssp_quad.ranges.size();
        if(in.ssp_quad.c_mat.size() != n_depth * n_range)
            throw Error("ssp_quad.c_mat 尺寸必须等于 ssp.points.size()*ssp_quad.ranges.size()");
    }

    if(in.boundaries.top_curve.r.size() != in.boundaries.top_curve.z.size() || in.boundaries.top_curve.r.size() < 2)
        throw Error("top_curve r/z 非法");
    if(in.boundaries.bottom_curve.r.size() != in.boundaries.bottom_curve.z.size() || in.boundaries.bottom_curve.r.size() < 2)
        throw Error("bottom_curve r/z 非法");
}

TL2DResult compute_tl_2d(const TL2DInput& in)
{
    validate(in);

    bhc::bhcInit init{};
    init.FileRoot = nullptr;

    // 回调：若用户没传，就默认打到 stderr（与 c_api 现行为一致）
    if(in.prt_callback) {
        init.prtCallback = +[](const char* m) { (void)m; }; // 占位，后面用捕获方式桥接（见下）
    }
    if(in.output_callback) {
        init.outputCallback = +[](const char* m) { (void)m; };
    }

    // 由于底层回调是 C 函数指针，不能直接传 std::function。
    // 这里采用“线程局部静态指针”做最小桥接（满足测试与单线程使用）。
    // 如需多实例并发，可升级为 map<thread_id, callbacks>。
    // 注意：MSVC 不允许在函数内本地类定义 static thread_local 成员。
    // 这里改用函数内的 thread_local 指针变量 + 无捕获静态函数桥接。
    static thread_local const TL2DInput* g_current_input = nullptr;

    struct CallbackBridge {
        static void prt(const char* m) {
            if(!m) return;
            if(g_current_input && g_current_input->prt_callback) g_current_input->prt_callback(m);
            else std::cerr << m;
        }
        static void out(const char* m) {
            if(!m) return;
            if(g_current_input && g_current_input->output_callback) g_current_input->output_callback(m);
            else std::cerr << m;
        }
    };

    init.prtCallback = &CallbackBridge::prt;
    init.outputCallback = &CallbackBridge::out;

    bhc::bhcParams<false> params;
    bhc::bhcOutputs<false, false> outputs;

    g_current_input = &in;
    g_progress_params = &params;
    g_progress_cache = 0;

    if(!bhc::setup_nofile<false, false>(init, params, outputs)) {
        g_current_input = nullptr;
        g_progress_params = nullptr;
        g_progress_cache = -1;
        throw Error("bhc::setup_nofile 失败");
    }

    // ---- 强制 no-file 模式：杜绝任何文件读取 ----
    // 与现 c_api/src/bhc_cpp_api.cpp 保持一致的防御性处理
    std::memset(params.Bdry->Top.hs.Opt, ' ', sizeof(params.Bdry->Top.hs.Opt));
    std::memset(params.Bdry->Bot.hs.Opt, ' ', sizeof(params.Bdry->Bot.hs.Opt));
    params.Bdry->Top.hs.Opt[0] = 'V';
    params.bdinfo->top.dirty = true;
    params.bdinfo->bot.dirty = true;

    // ---- 参数映射（全部来自 in 内存结构） ----
    set_title(params, in.title);
    set_freq0(params, in.freq_hz);
    set_runtype(params, in.run_type);
    set_beam(params, in.beam);
    set_positions(params, in.pos);
    set_angles(params, in.angles, in.n_rays, in.start_angle_deg, in.end_angle_deg);

    const bool is_quad_ssp = (in.options.find('Q') != std::string::npos);
    if(is_quad_ssp) {
        params.ssp->Type = 'Q';
        set_ssp_depth_grid(params, in.ssp);
        set_ssp_quad(params, in.ssp, in.ssp_quad);
    } else {
        params.ssp->Type = static_cast<char>(in.ssp.type);
        set_ssp_depth_grid(params, in.ssp);
    }

    set_boundaries_2d(params, in.boundaries);

    // ---- 运行 ----
    if(!bhc::echo<false>(params)) {
        bhc::finalize<false, false>(params, outputs);
        g_current_input = nullptr;
        g_progress_params = nullptr;
        g_progress_cache = -1;
        throw Error("bhc::echo 失败（请查看回调输出）");
    }

    if(!bhc::run<false, false>(params, outputs)) {
        bhc::finalize<false, false>(params, outputs);
        g_current_input = nullptr;
        g_progress_params = nullptr;
        g_progress_cache = -1;
        throw Error("bhc::run 失败（请查看回调输出）");
    }

    TL2DResult out;
    out.width  = params.Pos->NRr;
    out.height = params.Pos->NRz;
    out.tl_db.resize(static_cast<size_t>(out.width) * static_cast<size_t>(out.height));

    for(int ir = 0; ir < out.width; ++ir) {
        for(int iz = 0; iz < out.height; ++iz) {
            size_t idx = bhc::GetFieldAddr(0, 0, 0, 0, iz, ir, params.Pos);
            auto p     = outputs.uAllSources[idx];
            float amp  = std::sqrt(p.real() * p.real() + p.imag() * p.imag());
            out.tl_db[static_cast<size_t>(iz) * static_cast<size_t>(out.width) + static_cast<size_t>(ir)] =
                (amp > 0) ? -20.0f * std::log10(amp) : 200.0f;
        }
    }

    bhc::finalize<false, false>(params, outputs);
    g_current_input = nullptr;
    g_progress_params = nullptr;
    g_progress_cache = -1;
    return out;
}

} // namespace tedsbellhop

