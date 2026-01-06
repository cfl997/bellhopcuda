#include "TEDSbellhop/TEDSbellhop.hpp"

#include "bhc/bhc.hpp"   // 内部使用 bellhopcxx/bellhopcuda C++ API
#include "common.hpp"    // 与 c_api 保持一致（工程内已有）

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>

namespace tedsbellhop {

// =============================
// TL2DJob 内部实现
// =============================

struct TL2DJobImpl {
    std::atomic<int> progress{-1};
    std::atomic<bool> ready{false};

    // 用 mutex 保护 error / result / worker join
    mutable std::mutex m;
    std::string err;
    TL2DResult result;

    std::thread worker;

    // 供 worker 内部更新 progress：保存 params 指针
    bhc::bhcParams<false>* params_ptr = nullptr;

    ~TL2DJobImpl() {
        if(worker.joinable()) worker.join();
    }
};

TL2DJob::TL2DJob(TL2DJob&& o) noexcept = default;
TL2DJob& TL2DJob::operator=(TL2DJob&& o) noexcept = default;
TL2DJob::~TL2DJob() = default;

int TL2DJob::progress() const {
    if(!impl_) return -1;
    return impl_->progress.load(std::memory_order_relaxed);
}

bool TL2DJob::ready() const {
    if(!impl_) return true;
    return impl_->ready.load(std::memory_order_acquire);
}

std::string TL2DJob::error() const {
    if(!impl_) return {};
    std::lock_guard<std::mutex> lk(impl_->m);
    return impl_->err;
}

void TL2DJob::wait() {
    if(!impl_) return;
    if(impl_->worker.joinable()) impl_->worker.join();
}

TL2DResult TL2DJob::get() {
    wait();
    if(!impl_) return {};
    std::lock_guard<std::mutex> lk(impl_->m);
    if(!impl_->err.empty()) throw Error(impl_->err);
    return impl_->result;
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

    bhc::extsetup_sz<false>(params, static_cast<int32_t>(p.source_depths_m.size()));
    for(size_t i = 0; i < p.source_depths_m.size(); ++i) params.Pos->Sz[i] = p.source_depths_m[i];

    bhc::extsetup_rcvrranges<false>(params, static_cast<int32_t>(p.receiver_ranges.size()));
    for(size_t i = 0; i < p.receiver_ranges.size(); ++i) params.Pos->Rr[i] = p.receiver_ranges[i];

    bhc::extsetup_rcvrdepths<false>(params, static_cast<int32_t>(p.receiver_depths_m.size()));
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

    if(ssp_quad.c_mps.size() != static_cast<size_t>(n_depths) * static_cast<size_t>(n_ranges))
        throw Error("SSPQuad.c_mps 尺寸与深度/距离点数不匹配");

    bhc::extsetup_ssp_quad(params, n_depths, n_ranges);

    // 先设置单位标志，再写入 ranges 数据。
    // 说明：我们遇到过 ranges_in_km=true 时，bellhop 输出的 "Profile ranges (km):" 为空。
    // 这通常意味着 ranges 数组未被正确写入。为了稳健性，这里强制：
    // - 同时写 params.ssp->Seg.r 和 params.ssp->r（若存在）
    // - 并且写完后将 ssp 标记 dirty
    params.ssp->rangeInKm = ssp_quad.ranges_in_km;

    for(int i = 0; i < n_depths; ++i) {
        params.ssp->z[i] = ssp_base.points[static_cast<size_t>(i)].z_m;
    }

    for(int i = 0; i < n_ranges; ++i) {
        const auto rr = ssp_quad.ranges[static_cast<size_t>(i)];
        params.ssp->Seg.r[i] = rr;

    }

    // 布局与 c_api 当前实现一致：直接平铺复制
    for(size_t i = 0; i < ssp_quad.c_mps.size(); ++i) {
        params.ssp->cMat[i] = ssp_quad.c_mps[i];
    }

    params.ssp->dirty = true;
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

    //if(in.ssp.points.size() < 2) throw Error("ssp.points 至少需要2个点");
    for(size_t i = 1; i < in.ssp.points.size(); ++i) {
        if(!(in.ssp.points[i].z_m > in.ssp.points[i-1].z_m))
            throw Error("ssp.points 的 z_m 必须严格单调递增");
    }

    // Q-SSP 两种输入方式（优先级从高到低）：
    // 1) ssp_quad：SSP 网格独立于接收网格（ranges/depths/c_mps）
    // 2) 传统：options 包含 'Q' + ssp.type=Quad + ssp_quad

    const bool has_ssp_quad = !in.ssp_quad.c_mps.empty();
    const bool is_quad = has_ssp_quad || (in.options.find('Q') != std::string::npos);

    if(has_ssp_quad) {
        const size_t Nz = in.ssp_quad.depths_m.size();
        const size_t Nr = in.ssp_quad.ranges.size();
        if(Nz < 2 || Nr < 2) throw Error("使用 ssp_quad 时，ssp_quad.depths_m/ssp_quad.ranges 至少各需要2个点");
        if(in.ssp_quad.c_mps.size() != Nz * Nr) throw Error("ssp_quad.c_mps 尺寸必须等于 Nz*Nr");
    }

    if(!has_ssp_quad && is_quad) {
        if(in.ssp.type != SSPType::Quad) throw Error("options 包含 'Q' 时，ssp.type 必须为 SSPType::Quad");
        if(in.ssp_quad.ranges.size() < 2) throw Error("ssp_quad.ranges 至少需要2个点");
        const size_t n_depth = in.ssp.points.size();
        const size_t n_range = in.ssp_quad.ranges.size();
        if(in.ssp_quad.c_mps.size() != n_depth * n_range)
            throw Error("ssp_quad.c_mps 尺寸必须等于 ssp.points.size()*ssp_quad.ranges.size()");
    }

    if(in.boundaries.top_curve.r.size() != in.boundaries.top_curve.z.size() || in.boundaries.top_curve.r.size() < 2)
        throw Error("top_curve r/z 非法");
    if(in.boundaries.bottom_curve.r.size() != in.boundaries.bottom_curve.z.size() || in.boundaries.bottom_curve.r.size() < 2)
        throw Error("bottom_curve r/z 非法");
}

// 同步（阻塞）计算：内部复用异步实现，启动后等待 get()
TL2DResult compute_tl_2d(const TL2DInput& in)
{
    auto job = start_tl_2d(in);
    return job.get();
}


// 异步启动：新线程执行 bellhop，外部线程轮询 job.progress()
TL2DJob start_tl_2d(const TL2DInput& in)
{
    // 先在启动线程做输入校验（快速失败）
    validate(in);

    auto impl = std::make_shared<TL2DJobImpl>();
    impl->progress.store(0, std::memory_order_relaxed);

    // 注意：复制一份输入到 worker 线程，避免外部调用方修改 in 导致数据竞争
    // TL2DInput 内部含 vector/string/std::function，拷贝成本可接受（相对计算成本很小）。
    TL2DInput in_copy = in;

    impl->worker = std::thread([impl, in_copy = std::move(in_copy)]() mutable {
        // 为了支持从“外部线程”轮询进度：
        // - worker 线程里保存 params 指针到 impl->params_ptr
        // - 另起一个轻量循环在 worker 内更新 impl->progress（靠 bhc::get_percent_progress）

        try {
            bhc::bhcInit init{};
            init.FileRoot = nullptr;

            // ---- 回调桥接：每个 job 拥有自己的 input 副本 ----
            // bellhop 回调是 C 函数指针，这里用 thread_local 指针指向当前 job 的 input。
            // 注意：MSVC 不允许在函数内本地类声明 static data member，
            // 因此用“函数内 thread_local 变量 + lambda(无捕获)桥接”的方式。
            static thread_local const TL2DInput* g_current_input = nullptr;
            g_current_input = &in_copy;

            auto prt_cb = +[](const char* m) {
                if(!m) return;
                if(g_current_input && g_current_input->prt_callback) g_current_input->prt_callback(m);
                else std::cerr << m;
            };
            auto out_cb = +[](const char* m) {
                if(!m) return;
                if(g_current_input && g_current_input->output_callback) g_current_input->output_callback(m);
                else std::cerr << m;
            };

            init.prtCallback = prt_cb;
            init.outputCallback = out_cb;

            bhc::bhcParams<false> params;
            bhc::bhcOutputs<false, false> outputs;

            // 暴露 params 指针给进度查询
            impl->params_ptr = &params;

            if(!bhc::setup_nofile<false, false>(init, params, outputs)) {
                throw Error("bhc::setup_nofile 失败");
            }

            // ---- 强制 no-file 模式：杜绝任何文件读取 ----
            std::memset(params.Bdry->Top.hs.Opt, ' ', sizeof(params.Bdry->Top.hs.Opt));
            std::memset(params.Bdry->Bot.hs.Opt, ' ', sizeof(params.Bdry->Bot.hs.Opt));
            params.Bdry->Top.hs.Opt[0] = 'V';
            params.bdinfo->top.dirty = true;
            params.bdinfo->bot.dirty = true;

            // ---- 参数映射 ----
            set_title(params, in_copy.title);
            set_freq0(params, in_copy.freq_hz);
            set_runtype(params, in_copy.run_type);
            set_beam(params, in_copy.beam);
            set_positions(params, in_copy.pos);
            set_angles(params, in_copy.angles, in_copy.n_rays, in_copy.start_angle_deg, in_copy.end_angle_deg);

            // ---- SSP 映射 ----
            // 优先级（从高到低）：
            // 1) ssp_quad：SSP 网格独立于接收网格（你现在需要的模式）
            // 2) c_mps_grid：SSP 网格与接收网格完全一致
            // 3) 传统：按 options/ssp/ssp_quad

            const bool has_ssp_quad = !in_copy.ssp_quad.c_mps.empty();
            const bool want_quad_by_options = (in_copy.options.find('Q') != std::string::npos);

            if(has_ssp_quad) {
                const size_t Nz = in_copy.ssp_quad.depths_m.size();
                const size_t Nr = in_copy.ssp_quad.ranges.size();
                if(Nz < 2 || Nr < 2) throw Error("使用 ssp_quad 时，ssp_quad.depths_m/ssp_quad.ranges 至少各需要2个点");
                if(in_copy.ssp_quad.c_mps.size() != Nz * Nr)
                    throw Error("ssp_quad.c_mps 尺寸必须等于 Nz*Nr");

                params.ssp->Type = 'Q';

                SSP1D ssp_grid;
                ssp_grid.type = SSPType::Quad;
                ssp_grid.atten_unit = in_copy.ssp.atten_unit;
                ssp_grid.points.resize(Nz);
                for(size_t iz = 0; iz < Nz; ++iz) {
                    ssp_grid.points[iz].z_m = double(in_copy.ssp_quad.depths_m[iz]);
                    ssp_grid.points[iz].c_mps = double(in_copy.ssp_quad.c_mps[iz * Nr + 0]);
                }

                SSPQuad quad;
                quad.ranges_in_km = in_copy.ssp_quad.ranges_in_km;
                quad.ranges.resize(Nr);
                for(size_t ir = 0; ir < Nr; ++ir) quad.ranges[ir] = double(in_copy.ssp_quad.ranges[ir]);
                quad.c_mps.resize(Nz * Nr);
                for(size_t i = 0; i < Nz * Nr; ++i) quad.c_mps[i] = double(in_copy.ssp_quad.c_mps[i]);

                set_ssp_depth_grid(params, ssp_grid);
                set_ssp_quad(params, ssp_grid, quad);

            } else if(want_quad_by_options) {
                params.ssp->Type = 'Q';
                set_ssp_depth_grid(params, in_copy.ssp);
                set_ssp_quad(params, in_copy.ssp, in_copy.ssp_quad);
            } else {
                params.ssp->Type = static_cast<char>(in_copy.ssp.type);
                set_ssp_depth_grid(params, in_copy.ssp);
            }

            set_boundaries_2d(params, in_copy.boundaries);

            // ---- 在后台线程里更新进度（尽量低开销） ----
            // 说明：bhc::run 内部会更新进度；我们用另一个线程在 run 期间轮询。
            std::atomic<bool> run_done{false};
            std::thread progress_thread([&] {
                // 如果 get_percent_progress 在 preprocess 阶段也有意义，这里也能看到变化。
                while(!run_done.load(std::memory_order_acquire)) {
                    int p = bhc::get_percent_progress<false>(params);
                    impl->progress.store(p, std::memory_order_relaxed);
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                }
                // 结束时再刷一次
                int p = bhc::get_percent_progress<false>(params);
                impl->progress.store(p, std::memory_order_relaxed);
            });

            if(!bhc::echo<false>(params)) {
                run_done.store(true, std::memory_order_release);
                if(progress_thread.joinable()) progress_thread.join();
                bhc::finalize<false, false>(params, outputs);
                throw Error("bhc::echo 失败（请查看回调输出）");
            }

            if(!bhc::run<false, false>(params, outputs)) {
                run_done.store(true, std::memory_order_release);
                if(progress_thread.joinable()) progress_thread.join();
                bhc::finalize<false, false>(params, outputs);
                throw Error("bhc::run 失败（请查看回调输出）");
            }

            run_done.store(true, std::memory_order_release);
            if(progress_thread.joinable()) progress_thread.join();

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
            impl->progress.store(100, std::memory_order_relaxed);

            {
                std::lock_guard<std::mutex> lk(impl->m);
                impl->result = std::move(out);
            }

        } catch(const std::exception& e) {
            std::lock_guard<std::mutex> lk(impl->m);
            impl->err = e.what();
        }

        impl->params_ptr = nullptr;
        impl->ready.store(true, std::memory_order_release);
    });

    return TL2DJob{std::move(impl)};
}

} // namespace tedsbellhop

