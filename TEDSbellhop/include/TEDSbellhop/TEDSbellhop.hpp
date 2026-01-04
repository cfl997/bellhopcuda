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
struct SSPQuad {
    std::vector<double> ranges = {}; // m 或 km
    bool ranges_in_km = true;

    // 布局：range varies fastest
    // c_mat[iz * n_range + ir]
    std::vector<double> c_mat = {};
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
    SSP1D ssp;
    SSPQuad ssp_quad; // options 包含 'Q' 时使用

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
// API
// -----------------------------

TEDSBELLHOP_API void validate(const TL2DInput& in);

// -----------------------------
// 进度接口
// -----------------------------
// 说明：
// - 进度是“当前线程/当前任务”的百分比（0~100）。
// - 当没有正在运行的计算时，返回 -1。
// - 这是一个最小接口：为了做到像 background.cpp 那样轮询进度，
//   你需要用另一个线程调用 compute_tl_2d()，主线程里轮询 get_percent_progress()。
TEDSBELLHOP_API int get_percent_progress();

// 阻塞计算接口（与之前一致）
TEDSBELLHOP_API TL2DResult compute_tl_2d(const TL2DInput& in);

} // namespace tedsbellhop

