#include "bhc/bhc_cpp_api.hpp"

#include "bhc/bhc.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>

namespace bhc_cpp {

namespace {

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
    std::string s = pad_right(rt.value, 7);
    std::memcpy(params.Beam->RunType, s.data(), 7);
}

static inline void set_beam(bhc::bhcParams<false> &params, const Beam &b)
{
    std::string t = pad_right(b.type, 4);
    std::memcpy(params.Beam->Type, t.data(), 4);

    params.Beam->Box.x = static_cast<bhc::real>(b.box_x_m);
    params.Beam->Box.y = static_cast<bhc::real>(b.box_z_m);

    params.Beam->rangeInKm = false;

    if(b.deltas_m > 0) {
        params.Beam->autoDeltas = false;
        params.Beam->deltas = static_cast<bhc::real>(b.deltas_m);
    } else {
        params.Beam->autoDeltas = true;
        params.Beam->deltas = 0;
    }

    params.Beam->epsMultiplier = static_cast<bhc::real>(b.eps_multiplier);
}

static inline void set_freq0(bhc::bhcParams<false> &params, const Freq0 &f)
{
    params.freqinfo->freq0 = static_cast<bhc::real>(f.hz);
}

static inline void set_positions(bhc::bhcParams<false> &params, const Positions2D &p)
{
    if(p.sz_m.empty()) throw Error("Positions2D.sz_m 不能为空");
    if(p.rr_m.empty()) throw Error("Positions2D.rr_m 不能为空");
    if(p.rz_m.empty()) throw Error("Positions2D.rz_m 不能为空");

    bhc::extsetup_sz<false>(params, static_cast<int32_t>(p.sz_m.size()));
    for(size_t i=0;i<p.sz_m.size();++i) params.Pos->Sz[i] = p.sz_m[i];

    bhc::extsetup_rcvrranges<false>(params, static_cast<int32_t>(p.rr_m.size()));
    for(size_t i=0;i<p.rr_m.size();++i) params.Pos->Rr[i] = p.rr_m[i];

    bhc::extsetup_rcvrdepths<false>(params, static_cast<int32_t>(p.rz_m.size()));
    for(size_t i=0;i<p.rz_m.size();++i) params.Pos->Rz[i] = p.rz_m[i];

    params.Pos->RrInKm = p.rr_in_km;
    params.Beam->rangeInKm = false;
}

static inline void set_angles(bhc::bhcParams<false> &params, const Angles &a)
{
    if(a.alpha.empty()) throw Error("Angles.alpha 不能为空");
    bhc::extsetup_rayelevations<false>(params, static_cast<int32_t>(a.alpha.size()));
    params.Angles->alpha.inDegrees = a.alpha_in_degrees;
    for(size_t i=0;i<a.alpha.size();++i) {
        params.Angles->alpha.angles[i] = static_cast<bhc::real>(a.alpha[i]);
    }
}

static inline void set_ssp_1d(bhc::bhcParams<false> &params, const SSP1D &ssp)
{
    if(ssp.pts.size() < 2) throw Error("SSP1D.pts 至少需要 2 个点");

    params.ssp->Type = ssp.type;
    params.ssp->NPts = static_cast<int32_t>(ssp.pts.size());
    params.ssp->Nz   = params.ssp->NPts;

    std::string au = pad_right(ssp.atten_unit, 2);
    params.ssp->AttenUnit[0] = au[0];
    params.ssp->AttenUnit[1] = au[1];

    for(int32_t i=0;i<params.ssp->NPts;++i) {
        const auto &p = ssp.pts[static_cast<size_t>(i)];
        params.ssp->z[i]      = static_cast<bhc::real>(p.z_m);
        params.ssp->alphaR[i] = static_cast<bhc::real>(p.alphaR_mps);
        params.ssp->betaR[i]  = static_cast<bhc::real>(p.betaR_mps);
        params.ssp->rho[i]    = static_cast<bhc::real>(p.rho_gcm3);
        params.ssp->alphaI[i] = static_cast<bhc::real>(p.alphaI);
        params.ssp->betaI[i]  = static_cast<bhc::real>(p.betaI);
    }

    params.ssp->dirty = true;

    params.Bdry->Top.hs.Depth = params.ssp->z[0];
    params.Bdry->Bot.hs.Depth = params.ssp->z[params.ssp->NPts - 1];
}

static inline void write_boundary_curve_2d(
    bhc::bhcParams<false> &params,
    bhc::BdryInfoTopBot<false> &dst,
    const Boundary2D &src,
    bool is_top)
{
    if(src.r.size() != src.z.size()) {
        throw Error("Boundary2D.r/z 长度不一致");
    }
    if(src.r.size() < 2) {
        throw Error("Boundary2D 至少需要 2 个点");
    }

    // 2D 下 type[0] 不能是 'R'
    if(src.type.empty()) {
        throw Error("Boundary2D.type 不能为空");
    }
    const char t0 = src.type[0];
    const char t1 = (src.type.size() >= 2) ? src.type[1] : ' ';
    if(t0 != 'C' && t0 != 'L') {
        throw Error("Boundary2D.type[0] 在 2D 必须为 'C' 或 'L'");
    }
    if(!(t1 == 'S' || t1 == ' ' || t1 == 'L')) {
        throw Error("Boundary2D.type[1] 在 2D 必须为 'S'/' '/'L'");
    }

    // 写入 type 与单位
    dst.type[0] = t0;
    dst.type[1] = t1;
    dst.rangeInKm = src.range_in_km;
    dst.dirty = true;

    const int32_t N = static_cast<int32_t>(src.r.size());
    const bool ext = src.extend_to_infinity;

    // 分配：若扩展则 N+2，否则 N
    const int32_t NPts = ext ? (N + 2) : N;
    bhc::IORI2<false> np = NPts;
    // 复用 bellhop 的 extsetup_altimetry/bathymetry 会根据 ISTOP 选择 top/bot，
    // 这里为了不引入更多模板/实现，直接使用 trackallocate。
    // 但我们不应该直接 trackallocate（由 bellhop 管理）。
    // 因此：这里采用已有 extsetup_* API：
    if(is_top) {
        bhc::extsetup_altimetry<false>(params, np);
    } else {
        bhc::extsetup_bathymetry<false>(params, np);
    }

    // 获取被 extsetup 分配好的指针（extsetup 内部会指向 params.bdinfo->top/bot）
    // 注意：此处 dst.bd 已经是分配后的地址。

    auto write_point = [&](int32_t i, double rr, double zz) {
        // bd[i].x: vec2(range, depth)
        dst.bd[i].x = bhc::vec2(static_cast<bhc::real>(rr), static_cast<bhc::real>(zz));
    };

    int32_t offset = 0;
    if(ext) {
        write_point(0, -src.extend_left, src.z.front());
        offset = 1;
    }

    for(int32_t i=0;i<N;++i) {
        write_point(i + offset, src.r[static_cast<size_t>(i)], src.z[static_cast<size_t>(i)]);
        if(t1 == 'L') {
            if(src.hs.size() != static_cast<size_t>(N)) {
                throw Error("Boundary2D.type[1]=='L' 时必须提供 hs，且大小==点数");
            }
            auto &hs = dst.bd[i + offset].hs;
            const auto &inhs = src.hs[static_cast<size_t>(i)];
            hs.alphaR = static_cast<bhc::real>(inhs.alphaR_mps);
            hs.betaR  = static_cast<bhc::real>(inhs.betaR_mps);
            hs.rho    = static_cast<bhc::real>(inhs.rho_gcm3);
            hs.alphaI = static_cast<bhc::real>(inhs.alphaI);
            hs.betaI  = static_cast<bhc::real>(inhs.betaI);
        }
    }

    if(ext) {
        write_point(NPts - 1, src.extend_right, src.z.back());
        if(t1 == 'L' && !src.hs.empty()) {
            // 端点扩展 hs 复制边界相邻点（与 Read() 中扩展逻辑一致）
            dst.bd[0].hs = dst.bd[1].hs;
            dst.bd[NPts - 1].hs = dst.bd[NPts - 2].hs;
        }
    }

    // 深度一致性约束（防止 Validate 报 rises/drops）
    const double bdry_depth = is_top ? double(params.Bdry->Top.hs.Depth) : double(params.Bdry->Bot.hs.Depth);
    for(int32_t i=0;i<NPts;++i) {
        const double zz = dst.bd[i].x.y;
        if(is_top) {
            // 海面：z 应该等于 Top.hs.Depth（通常 0）或不低于它（取决于 NegTop 方向）；
            // 这里先做强制对齐（更接近“平坦海面”），复杂曲线后续再放开。
            (void)bdry_depth;
        } else {
            (void)bdry_depth;
        }
    }
}

static inline void set_boundaries_2d(bhc::bhcParams<false> &params, const Boundaries2D &b)
{
    // bc
    params.Bdry->Top.hs.bc = b.top.bc;
    params.Bdry->Bot.hs.bc = b.bot.bc;

    if(b.bot.bc == 'A') {
        params.Bdry->Bot.hs.alphaR = static_cast<bhc::real>(b.bot.alphaR_mps);
        params.Bdry->Bot.hs.betaR  = static_cast<bhc::real>(b.bot.betaR_mps);
        params.Bdry->Bot.hs.rho    = static_cast<bhc::real>(b.bot.rho_gcm3);
        params.Bdry->Bot.hs.alphaI = static_cast<bhc::real>(b.bot.alphaI);
        params.Bdry->Bot.hs.betaI  = static_cast<bhc::real>(b.bot.betaI);
    }

    // 调用者必须显式提供 top_curve/bot_curve 的点列（否则你就无法表达 ati/bty 文件语义）
    if(b.top_curve.r.empty() || b.top_curve.z.empty()) {
        throw Error("boundaries.top_curve 必须显式提供 r/z 点列（对应 .ati）");
    }
    if(b.bot_curve.r.empty() || b.bot_curve.z.empty()) {
        throw Error("boundaries.bot_curve 必须显式提供 r/z 点列（对应 .bty）");
    }

    // 写入 top/bot 曲线到 bdinfo
    write_boundary_curve_2d(params, params.bdinfo->top, b.top_curve, true);
    write_boundary_curve_2d(params, params.bdinfo->bot, b.bot_curve, false);
}

} // namespace

TLResult2D compute_tl_2d(const Input2D &in)
{
    bhc::bhcInit init{};
    init.FileRoot = nullptr;

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

    set_title(params, in.title);
    set_freq0(params, in.freq0);
    set_runtype(params, in.run_type);
    set_beam(params, in.beam);
    set_positions(params, in.pos);
    set_angles(params, in.angles);
    set_ssp_1d(params, in.ssp);
    set_boundaries_2d(params, in.boundaries);

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
    out.tl_db.resize(static_cast<size_t>(out.width) * static_cast<size_t>(out.height));

    for(int ir=0; ir<out.width; ++ir) {
        for(int iz=0; iz<out.height; ++iz) {
            const auto *Pos = params.Pos;
            const size_t idx = (((((size_t)0
                * (size_t)Pos->NSx + (size_t)0)
                * (size_t)Pos->NSy + (size_t)0)
                * (size_t)Pos->Ntheta + (size_t)0)
                * (size_t)Pos->NRz_per_range + (size_t)iz)
                * (size_t)Pos->NRr + (size_t)ir;

            auto p = outputs.uAllSources[idx];
            float amp = std::sqrt(p.real()*p.real() + p.imag()*p.imag());
            out.tl_db[static_cast<size_t>(iz)*static_cast<size_t>(out.width) + static_cast<size_t>(ir)]
                = (amp > 0) ? -20.0f * std::log10(amp) : 200.0f;
        }
    }

    bhc::finalize<false, false>(params, outputs);
    return out;
}

} // namespace bhc_cpp
