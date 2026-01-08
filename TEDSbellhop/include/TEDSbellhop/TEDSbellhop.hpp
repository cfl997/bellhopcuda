#pragma once

/**
 * @file TEDSbellhop.hpp
 * @brief TEDSbellhop：对 bellhopcxx / bellhopcuda 的“纯内存输入、C++17 友好”封装。
 *
 * ## 设计目标
 * - **只暴露 tedsbellhop 命名空间内的类型**：对外不直接暴露底层 bellhop 的 bhc:: 结构体/参数，避免 ABI/版本耦合。
 * - **不读取任何输入文件**：强制走 no-file 模式（所有环境、边界、SSP、射线、网格参数均由内存结构体提供）。
 * - **同步 + 异步两套调用方式**：既能阻塞式计算，也能启动后台任务并轮询进度。
 *
 * ## 使用方式概览
 * 1. 构造并填充 @ref tedsbellhop::TL2DInput
 * 2. （可选）调用 @ref tedsbellhop::validate 做参数校验（也可以让 compute/start 在内部校验）
 * 3. 调用：
 *    - 同步：@ref tedsbellhop::compute_tl_2d
 *    - 异步：@ref tedsbellhop::start_tl_2d -> @ref tedsbellhop::TL2DJob
 *
 * ## 单位与网格约定（非常重要）
 * - 深度统一使用 **米 (m)**。
 * - 距离（水平 range）既可用 m 也可用 km：
 *   - TL 输出网格使用 @ref tedsbellhop::Positions2D::receiver_ranges + @ref tedsbellhop::Positions2D::receiver_ranges_in_km
 *   - Q-SSP 的 range 网格使用 @ref tedsbellhop::SSPQuad::ranges + @ref tedsbellhop::SSPQuad::ranges_in_km
 * - Q-SSP 声速矩阵布局：@ref tedsbellhop::SSPQuad::c_mps 采用 **行优先（depth 行）**：
 *   `c_mps[iz * Nr + ir]`
 *
 * ## 与 bellhop 选项字符串（OPTIONS / RunType）关系
 * - bellhop 传统通过 7 字符 RunType + OPTIONS 控制算法分支。
 * - 本库提供：
 *   - “安全易用”的枚举组合（@ref tedsbellhop::RunType）
 *   - 以及用于完全复刻原行为的 raw 字符串（@ref tedsbellhop::RunType::raw + use_raw）。
 * - @ref tedsbellhop::TL2DInput::options 对应 env 的 OPTIONS1（如 "QVMT"），但库保证不读取任何文件。
 *
 * ## 线程安全与并行
 * - 多个 @ref tedsbellhop::TL2DJob 可并行启动（典型用于多剖面/多源深度批处理）。
 * - 具体并行度与底层实现（CPU/CUDA）及其内部资源管理有关。
 *
 * ## 说明
 * 本库内部仍会链接到底层 bellhop 静态库实现（bellhopcxxstatic / bellhopcudastatic）。
 */

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

/**
 * @brief 运行时错误类型。
 *
 * - 本库对外统一抛出/返回该错误类型（或其 std::runtime_error 基类语义）。
 * - 典型来源：参数校验失败、底层 bellhop 执行失败、CUDA/内存资源错误等。
 */
struct Error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

// -----------------------------
// 基础枚举：尽量贴近 bellhop 语义，但保持“可用且易懂”
// -----------------------------

/**
 * @brief 传输损失（TL）场的相干性选项。
 *
 * 对应 bellhop RunType 的第 1 个字符（TL option）。
 * - Coherent：相干叠加（保留相位），适合窄带/相干处理。
 * - SemiCoherent：半相干（常用于对干涉进行一定平均）。
 * - Incoherent：非相干叠加（功率叠加），常用于多路径能量统计。
 */
enum class TLCoherence : char { Coherent = 'C', SemiCoherent = 'S', Incoherent = 'I' };

/**
 * @brief 声线影响区/到达场幅度的近似模型（Influence / Amplitude model）。
 *
 * 对应 bellhop RunType 的第 2 个字符。
 *
 * 这些选项控制：如何从射线几何得到接收点的场/强度贡献（例如 Cerveny 形式、Gaussian/hat 函数等）。
 * 不同选项对数值稳定性、聚焦/散焦附近的行为会有显著影响。
 */
enum class Influence : char {
    /// Cerveny 振幅（笛卡尔坐标）
    CervenyCartesian        = 'C',
    /// Cerveny 振幅（ray-centered 坐标）
    CervenyRayCentered      = 'R',
    /// 简化 Gaussian 影响区
    SimpleGaussian          = 'S',
    /// 几何声学 Gaussian（ray-centered）
    GeomGaussianRayCentered = 'b',
    /// 几何声学 Gaussian（Cartesian）
    GeomGaussianCartesian   = 'B',
    /// 几何声学 hat 函数（ray-centered）
    GeomHatRayCentered      = 'g',
    /// 几何声学 hat 函数（Cartesian）
    GeomHatCartesian        = 'G'
};

/**
 * @brief 声源模型。
 *
 * 对应 bellhop RunType 的相关字符（不同版本 bellhop 可能位置略有差异）。
 * - Point：点源（常见）。
 * - Line：线源（用于 2D 模型中模拟无限长线源/柱面扩展）。
 */
enum class SourceModel : char { Point = 'R', Line = 'X' };

/**
 * @brief 接收点网格类型。
 *
 * - Rectilinear：规则二维网格（通常是 range×depth 的直角网格）。
 * - Irregular：不规则网格（由离散点列表给出）。
 */
enum class ReceiverGrid : char { Rectilinear = 'R', Irregular = 'I' };

/**
 * @brief 维度标志。
 *
 * - D2：标准 2D（r-z）传播。
 * - Nx2D：N×2D（多个方位/多剖面类用法，具体取决于底层 bellhop 版本）。
 * - D3：3D（如果底层支持）。
 */
enum class DimensionFlag : char { D2 = ' ', Nx2D = '2', D3 = '3' };

/**
 * @brief 声速剖面（SSP, Sound Speed Profile）类型。
 *
 * 对应 bellhop 的 SSP 选项。
 * - CLinear：c(z) 分段线性插值。
 * - N2Linear：n^2(z) 分段线性插值（某些情况下更稳定）。
 * - CubicSpline：三次样条。
 * - PCHIP：分段三次 Hermite（保形/单调性更好）。
 * - Analytic：解析型 SSP（由底层解析公式定义）。
 * - Quad：Q-SSP（2D、随距离变化的 SSP 网格，见 @ref tedsbellhop::SSPQuad）。
 */
enum class SSPType : char {
    CLinear     = 'C',
    N2Linear    = 'N',
    CubicSpline = 'S',
    PCHIP       = 'P',
    Analytic    = 'A',
    Quad        = 'Q'
};

/**
 * @brief 边界条件类型（海面/海底/半空间）。
 *
 * - Vacuum：真空边界（常用于自由表面近似：压力释放）。
 * - Rigid：刚性边界（法向速度为 0 的近似）。
 * - Acoustic：声学半空间（允许透射/吸收等，由 @ref tedsbellhop::Halfspace 参数描述）。
 */
enum class BoundaryCondition : char { Vacuum = 'V', Rigid = 'R', Acoustic = 'A' };

/**
 * @brief 边界曲线插值方式。
 *
 * - Curvilinear：曲线插值（更平滑，取决于底层实现）。
 * - Linear：折线线性插值（相邻节点线性连接）。
 */
enum class BoundaryInterp : char { Curvilinear = 'C', Linear = 'L' };

/**
 * @brief 边界曲线标志位。
 *
 * 目前仅保留 None（占位），用于与底层 bellhop 的 flag 字段语义对齐，
 * 未来如需支持更复杂的边界标志可在此扩展。
 */
enum class BoundaryFlag : char { None = ' ' };

// -----------------------------
// 对外输入结构
// -----------------------------

/**
 * @brief bellhop RunType 的封装。
 *
 * bellhop 传统使用一个固定长度（通常 7 字符）的字符串来控制：
 * - TL 场的相干性
 * - Influence / 幅度模型
 * - 源模型、接收网格、维度等
 *
 * 本库提供两种方式：
 * 1. **枚举字段组合**（推荐）：使用 tl/infl/source/grid/dim 等字段，由库内部生成对应 run-type。
 * 2. **完全手写 raw**：当你需要与某个现有 bellhop 输入文件保持逐字符一致时可用。
 *
 * @note
 * - 若 @ref use_raw = true，则优先使用 @ref raw（调用方需要自行保证长度与字符合法性）。
 * - 若 @ref use_raw = false，则由枚举字段组合生成。
 */
struct RunType {
    /// 是否直接使用 raw 字符串。
    bool use_raw = false;

    /// 原始 run-type 字符串（典型为 7 chars）。仅在 use_raw=true 时生效。
    std::string raw = "";

    TLCoherence tl = TLCoherence::Coherent;
    Influence infl = Influence::GeomHatCartesian;
    SourceModel source = SourceModel::Point;
    ReceiverGrid grid = ReceiverGrid::Rectilinear;
    DimensionFlag dim = DimensionFlag::D2;
};

/**
 * @brief 射线追踪/波束相关设置。
 *
 * 该结构体主要对应 bellhop 的 Beam 相关参数（尤其是计算域与步长）。
 *
 * @note
 * - box_x_m / box_z_m 通常决定射线在数值积分时允许传播的最大范围。
 * - deltas_m <= 0 时，底层通常会选择一个“自动步长”。
 */
struct Beam {
    /// bellhop 中 Beam.Type 通常为 4 字符；本封装保留字符串以便兼容（常用 'G'/'B' 等）。
    std::string type = "G";

    /// 最大水平距离（单位 m）。
    double box_x_m = 10000.0;

    /// 最大深度（单位 m）。
    double box_z_m = 200.0;

    /// 射线积分步长（单位 m）。<=0 表示由底层自动选择。
    double deltas_m = -1.0;

    /// 某些底层实现会用该系数缩放 eps/容差（用于数值稳定或性能调参）。
    double eps_multiplier = 1.0;
};

/**
 * @brief 2D 几何位置输入（声源/接收点）。
 *
 * - 声源支持多个深度（可用于一次性计算多源深度的 TL；具体是否合并取决于底层实现）。
 * - 接收点网格通常是 range×depth 的直角网格：
 *   - receiver_ranges：水平距离采样点
 *   - receiver_depths_m：深度采样点
 *
 * @note
 * - 深度单位固定为 m。
 * - range 单位由 receiver_ranges_in_km 指定。
 */
struct Positions2D {
    /// 声源深度列表（单位 m）。
    std::vector<float> source_depths_m = {50.0f};

    /// 接收点水平距离列表（单位 m 或 km，见 receiver_ranges_in_km）。
    std::vector<float> receiver_ranges = {};

    /// 接收点深度列表（单位 m）。
    std::vector<float> receiver_depths_m = {};

    /// receiver_ranges 是否以 km 为单位。false 表示 m。
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
    // 距离节点（长度 range_count），单位由 ranges_in_km 决定
    // - 在 TL2DResult 中对应 width（水平采样点数）
    // - 在 c_mps 中对应“列数”（range 为内层连续维度）
    std::vector<float> ranges = {}; // m 或 km
    bool ranges_in_km = false;

    // 深度节点（长度 depth_count），单位 m
    // - 在 TL2DResult 中对应 height（深度采样点数）
    // - 在 c_mps 中对应“行数”（depth 为外层维度）
    std::vector<float> depths_m = {};

    // 声速矩阵（m/s），尺寸 Nz*Nr
    // 布局约定（row-major，depth 行优先；range 为内层连续维度）：
    //   - 令 width  = range_count  = ranges.size()  （列数）
    //   - 令 height = depth_count  = depths_m.size()（行数）
    //   - 则访问方式：c_mps[depth_index * width + range_index]
    std::vector<float> c_mps = {};
};

/**
 * @brief 声学半空间（Half-space）参数。
 *
 * 用于描述上边界/下边界的“介质属性”和边界条件。
 * 在 bellhop 语义中，边界（surface/bottom）通常由两部分组成：
 * 1) **边界几何**：由 @ref tedsbellhop::BoundaryCurve2D 给出 (r,z) 曲线
 * 2) **边界介质/条件**：由 @ref tedsbellhop::Halfspace 给出（刚性/真空/声学半空间等）
 *
 * @note
 * - 当 bc = @ref tedsbellhop::BoundaryCondition::Rigid 或 Vacuum 时，下方声学参数通常不会被使用。
 * - 当 bc = @ref tedsbellhop::BoundaryCondition::Acoustic 时，这些参数用于计算反射/透射、吸收等效应。
 *
 * 字段含义（尽量与 bellhop 命名保持一致）：
 * - alphaR_mps：P 波（纵波）**相速度实部**，单位 m/s
 * - betaR_mps ：S 波（横波）相速度实部，单位 m/s（流体中通常为 0）
 * - rho_gcm3  ：密度，单位 g/cm^3
 * - alphaI / betaI：速度虚部或与衰减相关的参数（具体解释依赖 bellhop 版本/单位系统）
 * - opt：半空间的附加选项字符串（最多 6 字符，原样传递给底层）
 * - sigma：界面粗糙度/统计参数（若底层启用对应模型时使用）
 */
struct Halfspace {
    /// 边界条件类型：Vacuum/Rigid/Acoustic。
    BoundaryCondition bc = BoundaryCondition::Rigid;

    /// 纵波速度实部 (m/s)。
    double alphaR_mps = 1600.0;

    /// 横波速度实部 (m/s)。流体介质可置 0。
    double betaR_mps  = 0.0;

    /// 密度 (g/cm^3)。
    double rho_gcm3   = 1.8;

    /// 纵波虚部/吸收参数（与 bellhop 的衰减模型/单位相关）。
    double alphaI     = 0.0;

    /// 横波虚部/吸收参数（与 bellhop 的衰减模型/单位相关）。
    double betaI      = 0.0;

    /// 额外选项（最多 6 字符）。仅在底层模型支持时生效。
    std::string opt = "";

    /// 粗糙度/散射相关参数（仅在底层启用相应模型时生效）。
    double sigma = 0.0;
};

/**
 * @brief 2D 边界曲线（surface 或 bottom）的几何描述。
 *
 * 使用一组 (r,z) 节点定义边界随距离变化的深度。
 *
 * - r：水平距离节点（单位由 range_in_km 指定）
 * - z：对应深度节点（单位 m）
 *
 * 常见用法：
 * - 海面：通常为 z=0 的水平线（也可以给出随距离起伏的海面）。
 * - 海底：给出随距离变化的地形曲线。
 *
 * @note
 * - r 与 z 的长度必须一致，且通常要求 r 单调递增。
 * - interp 控制节点间如何插值（线性/曲线）。
 * - extend_to_infinity=true 时，底层可将曲线在左右端点之外“延拓”，避免射线跑出定义域。
 *   extend_left/extend_right 为延拓距离（单位与 r 一致；默认给很大值相当于“无限延拓”）。
 */
struct BoundaryCurve2D {
    /// 插值方式：Linear/Curvilinear。
    BoundaryInterp interp = BoundaryInterp::Linear;

    /// 曲线标志位（目前仅 None，占位）。
    BoundaryFlag flag = BoundaryFlag::None;

    /// r 的单位是否为 km。false 表示 m。
    bool range_in_km = false;

    /// 水平距离节点（m 或 km）。
    std::vector<double> r = {};

    /// 深度节点（m）。
    std::vector<double> z = {};

    /// 是否在左右端点外进行延拓。
    bool extend_to_infinity = true;

    /// 左侧延拓长度（单位同 r，默认极大）。
    double extend_left = 1e9;

    /// 右侧延拓长度（单位同 r，默认极大）。
    double extend_right = 1e9;
};

/**
 * @brief 2D 场景的上下边界组合。
 *
 * - top / bottom：边界介质与边界条件（@ref tedsbellhop::Halfspace）
 * - top_curve / bottom_curve：边界几何曲线（@ref tedsbellhop::BoundaryCurve2D）
 *
 * @note
 * 常见设置：
 * - top.bc = Vacuum，top_curve 为 z=0 的水平线
 * - bottom.bc = Rigid 或 Acoustic，bottom_curve 为海底地形
 */
struct Boundaries2D {
    Halfspace top;
    Halfspace bottom;
    BoundaryCurve2D top_curve;
    BoundaryCurve2D bottom_curve;
};

/**
 * @brief 2D 传输损失（TL）计算的完整输入。
 *
 * 该结构体等价于 bellhop 的一整套环境输入（env、bty、ssp、射线、接收网格等），
 * 但完全以内存方式提供，不涉及任何文件读写。
 *
 * ### 必填/常用字段
 * - title：标题（用于日志/调试标识）
 * - freq_hz：频率（Hz）
 * - pos.receiver_ranges / pos.receiver_depths_m：TL 输出网格
 * - pos.source_depths_m：声源深度（m）
 * - boundaries：上下边界（介质 + 几何）
 * - ssp 或 ssp_quad：声速剖面（1D 或 2D Q-SSP）
 *
 * ### OPTIONS 与 Q-SSP 选择
 * - options 对应 bellhop env 的 OPTIONS1（例如 "QVMT"）。
 * - 对于 2D Q-SSP：
 *   - 推荐：填充 @ref ssp_quad.depths_m / @ref ssp_quad.ranges / @ref ssp_quad.c_mps
 *   - 同时在 options 中包含 'Q'
 *   - 某些实现中即使不写 'Q'，库也可能在检测到 c_mps 非空时强制走 Q-SSP（具体由实现决定）。
 *
 * ### 射线角度生成规则
 * - 若 angles.alpha 非空：直接使用 angles.alpha（单位由 angles.alpha_in_degrees 指定）
 * - 若 angles.alpha 为空：使用 n_rays + start_angle_deg/end_angle_deg 生成等间隔角度
 *
 * ### 日志回调
 * - prt_callback / output_callback 可接收底层输出文本（用于嵌入式 UI/IDE 日志面板）。
 * - 回调参数为 `const char*`，生命周期仅在回调调用期间有效（调用方不要缓存指针）。
 */
struct TL2DInput {
    /// 标题/任务名（用于日志与输出标识）。
    std::string title = "";

    /// 频率（Hz）。
    double freq_hz = 1000.0;

    /// 对应 bellhop env 的 OPTIONS1（例如 "QVMT"）。本库保证不读取任何文件。
    std::string options = "";

    /// 当 angles.alpha 为空时，用该参数生成射线角数量。
    int n_rays = 0;

    /// 当 angles.alpha 为空时，生成射线角的起始角（度）。
    double start_angle_deg = -90.0;

    /// 当 angles.alpha 为空时，生成射线角的终止角（度）。
    double end_angle_deg = 90.0;

    /// RunType（推荐使用枚举字段；必要时可 use_raw）。
    RunType run_type;

    /// Beam/射线积分与计算域参数。
    Beam beam;

    /// 声源/接收点网格。
    Positions2D pos;

    /// 射线初始角设置。
    Angles angles;

    /// 1D 声速剖面（默认路径）。
    SSP1D ssp;

    /// 2D Q-SSP（range-dependent SSP）。
    SSPQuad ssp_quad;

    /// 上下边界（海面/海底）的介质与几何。
    Boundaries2D boundaries;

    /// 可选日志回调：对应底层“打印/日志”通道。
    std::function<void(const char*)> prt_callback;

    /// 可选日志回调：对应底层“输出/进度/结果摘要”通道。
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
