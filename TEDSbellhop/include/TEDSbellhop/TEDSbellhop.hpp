#pragma once

// TEDSbellhop - C++17 friendly wrapper library for bellhopcxx/bellhopcuda
// 目标：
//  - 对外只暴露 TEDSbellhop 自己的结构体（不暴露任何 bhc:: 参数/结构体）
//  - 不读取任何输入文件（强制 no-file 模式）
//  - 通过内存结构体提供全部参数并计算
//
// 说明：本库内部仍会链接到底层 bellhop 静态库实现（bellhopcxxstatic/bellhopcudastatic）。

#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_WIN32)
    #if defined(TEDSBELLHOP_EXPORTS)
        #define TEDSBELLHOP_API __declspec(dllexport)
    #else
        #define TEDSBELLHOP_API __declspec(dllimport)
    #endif
#else
    #define TEDSBELLHOP_API __attribute__((visibility("default")))
#endif

namespace tedsbellhop {

struct Error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

// -----------------------------
// 基础枚举：尽量贴近 bellhop 语义，但保持“可用且易懂”
// -----------------------------

enum class TLCoherence : char { Coherent = 'C', SemiCoherent = 'S', Incoherent = 'I' };

enum class Influence : char {
    CervenyCartesian        = 'C',
    CervenyRayCentered      = 'R',
    SimpleGaussian          = 'S',
    GeomGaussianRayCentered = 'b',
    GeomGaussianCartesian   = 'B',
    GeomHatRayCentered      = 'g',
    GeomHatCartesian        = 'G'
};

enum class SourceModel : char { Point = 'R', Line = 'X' };

enum class ReceiverGrid : char { Rectilinear = 'R', Irregular = 'I' };

enum class DimensionFlag : char { D2 = ' ', Nx2D = '2', D3 = '3' };

enum class SSPType : char {
    CLinear     = 'C',
    N2Linear    = 'N',
    CubicSpline = 'S',
    PCHIP       = 'P',
    Analytic    = 'A',
    Quad        = 'Q'
};

enum class BoundaryCondition : char { Vacuum = 'V', Rigid = 'R', Acoustic = 'A' };

enum class BoundaryInterp : char { Curvilinear = 'C', Linear = 'L' };

enum class BoundaryFlag : char { None = ' ' };

// -----------------------------
// 对外输入结构
// -----------------------------

struct RunType {
    // 提供“安全枚举组合”，也保留 raw 以便完全复刻 bellhop 7字符行为
    bool use_raw = false;
    std::string raw = ""; // 7 chars

    TLCoherence tl = TLCoherence::Coherent;
    Influence infl = Influence::GeomHatCartesian;
    SourceModel source = SourceModel::Point;
    ReceiverGrid grid = ReceiverGrid::Rectilinear;
    DimensionFlag dim = DimensionFlag::D2;
};

struct Beam {
    // bellhop 里 Beam.Type 为 4 字符；这里默认给 "G   "（帽函数相关）
    std::string type = "G";
    double box_x_m = 10000.0; // 最大水平距离
    double box_z_m = 200.0;   // 最大深度
    double deltas_m = -1.0;   // <=0 自动
    double eps_multiplier = 1.0;
};

struct Positions2D {
    std::vector<float> source_depths_m = {50.0f};
    std::vector<float> receiver_ranges = {}; // m 或 km
    std::vector<float> receiver_depths_m = {};
    bool receiver_ranges_in_km = false;
};

struct Angles {
    // 若为空，将用 n_rays/start/end 自动生成
    std::vector<double> alpha = {};
    bool alpha_in_degrees = true;
};

struct SSPPoint {
    double z_m = 0.0;
    double c_mps = 1500.0;
    double cs_mps = 0.0;
    double rho_gcm3 = 1.0;
    double atten_c = 0.0;
    double atten_s = 0.0;
};

struct SSP1D {
    SSPType type = SSPType::CLinear;
    std::string atten_unit = "W ";
    std::vector<SSPPoint> points = {}; // 至少2个，z 单调递增
};

// Type 'Q' 的 2D range-dependent SSP
// 说明：
// - bellhop Q-SSP 的网格由 (depths_m, ranges) 定义。
// - TL 的输出网格（Positions2D.receiver_ranges/receiver_depths_m）可以与 Q-SSP 网格不同。
struct SSPQuad {
    // 距离节点（长度 Nr），单位由 ranges_in_km 决定
    std::vector<float> ranges = {}; // m 或 km
    bool ranges_in_km = false;

    // 深度节点（长度 Nz），单位 m
    std::vector<float> depths_m = {};

    // 声速矩阵（m/s），尺寸 Nz*Nr
    // 布局约定（range 优先）：c_mps[iz * Nr + ir]
    std::vector<float> c_mps = {};
};

struct Halfspace {
    BoundaryCondition bc = BoundaryCondition::Rigid;

    double alphaR_mps = 1600.0;
    double betaR_mps  = 0.0;
    double rho_gcm3   = 1.8;
    double alphaI     = 0.0;
    double betaI      = 0.0;

    std::string opt = ""; // up to 6 chars
    double sigma = 0.0;
};

struct BoundaryCurve2D {
    BoundaryInterp interp = BoundaryInterp::Linear;
    BoundaryFlag flag = BoundaryFlag::None;

    bool range_in_km = false;
    std::vector<double> r = {};
    std::vector<double> z = {};

    bool extend_to_infinity = true;
    double extend_left = 1e9;
    double extend_right = 1e9;
};

struct Boundaries2D {
    Halfspace top;
    Halfspace bottom;
    BoundaryCurve2D top_curve;
    BoundaryCurve2D bottom_curve;
};

struct TL2DInput {
    std::string title = "";
    double freq_hz = 1000.0;

    // 对应 env 的 OPTIONS1（例如 "QVMT"），但本库保证不读取任何文件
    std::string options = "";

    // 当 angles.alpha 为空时，用下面参数生成等间隔射线角
    int n_rays = 0;
    double start_angle_deg = -90.0;
    double end_angle_deg = 90.0;

    RunType run_type;
    Beam beam;
    Positions2D pos;
    Angles angles;

    // SSP
    // 说明：
    // - 默认走 1D SSP：ssp
    // - 若要使用 2D Q-SSP：
    //   1) 填写 ssp_quad.depths_m / ssp_quad.ranges / ssp_quad.c_mps
    //   2) options 里包含 'Q'（或直接留空，库内部在检测到 ssp_quad.c_mps 非空时会强制走 Q-SSP）
    SSP1D ssp;
    SSPQuad ssp_quad;

    Boundaries2D boundaries;

    // 可选日志回调（不强制）
    std::function<void(const char*)> prt_callback;
    std::function<void(const char*)> output_callback;
};

struct TL2DResult {
    int width = 0;  // ranges count
    int height = 0; // depths count
    std::vector<float> tl_db; // tl_db[iz*width+ir]

    float at(int ir, int iz) const {
        return tl_db.at(static_cast<size_t>(iz) * static_cast<size_t>(width) + static_cast<size_t>(ir));
    }
};

// -----------------------------
// 异步任务（Job）接口
// -----------------------------

struct TL2DJobImpl;

class TEDSBELLHOP_API TL2DJob {
public:
    TL2DJob() = default;

    // 可移动，不可拷贝
    TL2DJob(TL2DJob&&) noexcept;
    TL2DJob& operator=(TL2DJob&&) noexcept;
    TL2DJob(const TL2DJob&) = delete;
    TL2DJob& operator=(const TL2DJob&) = delete;

    ~TL2DJob();

    // 当前进度：0~100；未开始/不可用返回 -1
    int progress() const;

    // 是否完成（成功或失败都算完成）
    bool ready() const;

    // 若失败，返回错误信息；未失败返回空串
    std::string error() const;

    // 等待完成并取结果；如果任务失败会抛 tedsbellhop::Error
    TL2DResult get();

    // 主动等待完成（不取结果）
    void wait();

    // 是否持有有效任务
    explicit operator bool() const noexcept { return static_cast<bool>(impl_); }

private:
    friend TEDSBELLHOP_API TL2DJob start_tl_2d(const TL2DInput& in);
    explicit TL2DJob(std::shared_ptr<TL2DJobImpl> impl) : impl_(std::move(impl)) {}

    std::shared_ptr<TL2DJobImpl> impl_;
};

// -----------------------------
// API
// -----------------------------

TEDSBELLHOP_API void validate(const TL2DInput& in);

// 同步（阻塞）接口
TEDSBELLHOP_API TL2DResult compute_tl_2d(const TL2DInput& in);

// 异步接口：启动计算并返回 job
TEDSBELLHOP_API TL2DJob start_tl_2d(const TL2DInput& in);

} // namespace tedsbellhop
