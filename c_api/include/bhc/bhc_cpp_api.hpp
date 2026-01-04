#pragma once

// 说明：这是 bellhopcxx/bellhopcuda 的 C++17 友好封装接口。
// 目标：完全不读任何输入文件，外部以 C++ 结构体形式提供全部参数。

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------
// DLL 导出/导入宏
// -----------------------------
#if defined(_WIN32)
    #if defined(BHC_CPP_API_EXPORTS)
        #define BHC_CPP_API __declspec(dllexport)
    #else
        #define BHC_CPP_API __declspec(dllimport)
    #endif
#else
    #define BHC_CPP_API __attribute__((visibility("default")))
#endif

namespace bhc_cpp {

// -----------------------------
// 基础工具和枚举类型
// -----------------------------

/** @brief API 抛出的标准异常类型 */
struct Error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/** @enum SSPType 声速剖面插值类型 (对应 env 文件 SSP 段的选项) */
enum class SSPType : char {
    CLinear     = 'C', /**< C-线性插值 */
    N2Linear    = 'N', /**< N²-线性插值 */
    CubicSpline = 'S', /**< 三次样条插值 */
    PCHIP       = 'P', /**< PCHIP (分段三次Hermite) 插值 */
    Analytic    = 'A', /**< 解析式声速剖面 */
    Quad        = 'Q'  /**< 2D Range-Dependent SSP (Quadrilateral) */
};

/** @enum BoundaryCondition 海面/海底边界条件类型 (对应 env 文件 TopOpt/BotOpt 的 BC 选项) */
enum class BoundaryCondition : char {
    Vacuum   = 'V', /**< 真空 (完全反射) */
    Rigid    = 'R', /**< 刚性 (声压梯度为0) */
    Acoustic = 'A'  /**< 声学半空间 (需要提供底层介质参数) */
};

/** @enum SBPFlag 源波束方向图选项 (对应 env 文件 SBP 段的 SBPFLG) */
enum class SBPFlag : char {
    None = 'N', /**< 无方向图 (全向声源) */
    File = 'F'  /**< 使用内存中的 vector 数据 (等价于从 .sbp 文件读取) */
};

/** @enum TLCoherence TL 相干性 (对应 RunType[0]) */
enum class TLCoherence : char {
    Coherent     = 'C', /**< 相干 TL */
    SemiCoherent = 'S', /**< 半相干 TL */
    Incoherent   = 'I'  /**< 非相干 TL */
};

/** @enum Influence 影响函数/波束类型 (对应 RunType[1]) */
enum class Influence : char {
    CervenyCartesian        = 'C', /**< Cerveny Cartesian beams */
    CervenyRayCentered      = 'R', /**< Cerveny ray centered beams */
    SimpleGaussian          = 'S', /**< Simple gaussian beams */
    GeomGaussianRayCentered = 'b', /**< Geometric gaussian beams in ray-centered coordinates */
    GeomGaussianCartesian   = 'B', /**< Geometric gaussian beams in Cartesian coordinates */
    GeomHatRayCentered      = 'g', /**< Geometric hat beams in ray-centered coordinates */
    GeomHatCartesian        = 'G'  /**< Geometric hat beams in Cartesian coordinates */
};

/** @enum SourceModel 源模型 (对应 RunType[3]) */
enum class SourceModel : char {
    Point = 'R', /**< 点源 (柱坐标) */
    Line  = 'X'  /**< 线源 (笛卡尔坐标) */
};

/** @enum ReceiverGrid 接收器网格 (对应 RunType[4]) */
enum class ReceiverGrid : char {
    Rectilinear = 'R', /**< 规则网格 */
    Irregular   = 'I'  /**< 不规则网格 */
};

/** @enum DimensionFlag 维度标志 (对应 RunType[5]) */
enum class DimensionFlag : char {
    D2   = ' ', /**< 2D 计算 */
    Nx2D = '2', /**< N x 2D 计算 */
    D3   = '3'  /**< 3D 计算 */
};

/** @struct RunTypeEnum 以枚举形式表达 RunType 7 字符 */
struct RunTypeEnum {
    TLCoherence tl     = TLCoherence::Coherent;
    Influence infl     = Influence::GeomHatCartesian;
    SourceModel source = SourceModel::Point;
    ReceiverGrid grid  = ReceiverGrid::Rectilinear;
    DimensionFlag dim  = DimensionFlag::D2;
};

// -----------------------------
// 与 "文件语义" 等价的输入结构
// -----------------------------

/** @struct Title 对应 env 文件中的标题行 */
struct Title {
    std::string text = ""; /**< 运行标题，将被截断或填充到80字符 */
};

/** @struct Freq0 对应 env 文件中的频率设置 */
struct Freq0 {
    double hz = 1000.0; /**< 中心频率 (Hz) */
};

/** @struct RunType 对应 env 文件中的 RunType 字符串 */
struct RunType {
    bool use_raw    = false; /**< 为 true 时使用下面的 raw 字符串，否则使用 e 枚举组合 */
    std::string raw = "";    /**< 原始7字符 RunType 字符串，用于完全复刻 env 文件行为 */
    RunTypeEnum e{};         /**< 类型安全的枚举组合，默认生效 */
};

/** @struct Beam 对应 env 文件中的 Beam/Box 设置 */
struct Beam {
    std::string type      = "G";     /**< 影响函数类型 (4字符)，与 RunType[1] 配合使用 */
    double box_x_m        = 10000.0; /**< 射线追踪的最大水平距离 (m) */
    double box_z_m        = 200.0;   /**< 射线追踪的最大深度 (m) */
    double deltas_m       = -1.0;    /**< 射线步长 (m)，<=0 表示自动选择；对应 env 里 RAY-STEP */
    double eps_multiplier = 1.0;     /**< 用于调整波束水平宽度的乘子 */
};

/** @struct Positions2D 对应 env 文件中的源/接收器位置 (2D) */
struct Positions2D {
    std::vector<float> sz_m = {50.0f}; /**< 源深度数组 (m) */
    std::vector<float> rr_m = {};      /**< 接收器距离数组 (m 或 km，由 rr_in_km 控制) */
    std::vector<float> rz_m = {};      /**< 接收器深度数组 (m) */
    bool rr_in_km           = false;   /**< 若为 true，则 rr_m 中的单位为公里 */
};

/** @struct Angles 对应 env 文件中的射线出射角 */
struct Angles {
    std::vector<double> alpha = {};   /**< 射线出射仰角数组 */
    bool alpha_in_degrees     = true; /**< 若为 true，则 alpha 中的单位为度，否则为弧度 */
};

/** @struct SSP1DPoint 对应 env 文件中 SSP 段的单行数据 */
struct SSP1DPoint {
    double z_m        = 0.0;    /**< 深度 (m) */
    double alphaR_mps = 1500.0; /**< 压缩波声速 (m/s) */
    double betaR_mps  = 0.0;    /**< 剪切波声速 (m/s)，水中为0 */
    double rho_gcm3   = 1.0;    /**< 密度 (g/cm³) */
    double alphaI     = 0.0;    /**< 压缩波衰减 (单位由 SSP1D.atten_unit 决定) */
    double betaI      = 0.0;    /**< 剪切波衰减 (单位由 SSP1D.atten_unit 决定) */
};

/** @struct SSP1D 对应 env 文件中的主 SSP 段 */
struct SSP1D {
    SSPType type           = SSPType::CLinear; /**< 插值类型 */
    std::string atten_unit = "W ";    /**< 衰减单位 (2字符)，如 "W "(dB/wavelength), "F "(Thorp)等 */
    std::vector<SSP1DPoint> pts = {}; /**< 声速剖面点，必须至少2个，且深度 z_m 须单调递增 */
};

/** @struct SSPQuad 2D Range-Dependent SSP (Quadrilateral, Type 'Q') */
struct SSPQuad {
    std::vector<double> ranges = {}; /**< 各个声速剖面所在的距离 (单位由 ranges_in_km 控制) */
    bool ranges_in_km = true;        /**< ranges 向量的单位是否为公里 */
    /**
     * @brief 声速矩阵 (m/s)，平铺存储。
     * 尺寸: ssp1d.pts.size() * ranges.size()
     * 布局: range 优先（range varies fastest）。
     * 即 ssp_matrix[i_depth * ranges.size() + i_range]
     */
    std::vector<double> ssp_matrix = {};
};

/** @enum BoundaryInterp 2D 边界曲线插值方式（对应 .ati/.bty 的首字符） */
enum class BoundaryInterp : char {
    Curvilinear = 'C', /**< Curvilinear / spline-like */
    Linear      = 'L'  /**< Piecewise linear */
};

/**
 * @enum BoundaryFlag 2D 边界曲线的第二字符标志位。
 * 说明：bellhop 的第二字符在不同版本/模式下语义可能不同；
 * 本封装默认不强行解释，只透传到内部。
 */
enum class BoundaryFlag : char {
    None = ' ' /**< 无标志 */
};

/** @struct Boundary2D 对应 .ati/.bty 文件的2D边界曲线 */
struct Boundary2D {
    BoundaryInterp interp = BoundaryInterp::Linear;
    BoundaryFlag flag     = BoundaryFlag::None;

    bool range_in_km      = false;   /**< 若为 true，则 r 中的单位为公里 */
    std::vector<double> r = {};      /**< 距离点，必须单调递增，至少2个点 */
    std::vector<double> z = {};      /**< 深度点，与 r 等长 */
    std::vector<SSP1DPoint> hs = {}; /**< 预留：当某些模式需要为每个点提供半空间参数 */

    bool extend_to_infinity = true;  /**< 是否自动将边界扩展到 +/- 无穷远 */
    double extend_left      = 1e9;   /**< 左侧扩展距离 (m 或 km，取决于 range_in_km) */
    double extend_right     = 1e9;   /**< 右侧扩展距离 (m 或 km，取决于 range_in_km) */
};

/** @struct Halfspace 对应 env 文件中的声学半空间参数 */
struct Halfspace {
    BoundaryCondition bc = BoundaryCondition::Rigid; /**< 边界条件类型 */

    // 对应 HSInfo
    double alphaR_mps = 1600.0; /**< 半空间压缩波声速 (m/s) */
    double betaR_mps  = 0.0;    /**< 半空间剪切波声速 (m/s) */
    double rho_gcm3   = 1.8;    /**< 半空间密度 (g/cm³) */
    double alphaI     = 0.0;    /**< 半空间压缩波衰减 */
    double betaI      = 0.0;    /**< 半空间剪切波衰减 */

    // 对应 HSInfo::Opt[6] 与 HSExtra::Sigma
    std::string opt = ""; /**< TopOpt/BotOpt 附加选项（0~6 字符，将右侧填空到 6 字符） */
    double sigma    = 0.0; /**< 粗糙度/散射相关参数（bellhop 读入并 echo，通常不参与物理计算） */
};

/** @struct Boundaries2D 对应 env 文件中的 TopOpt/BotOpt 及 .ati/.bty */
struct Boundaries2D {
    Halfspace top;        /**< 海面半空间参数 */
    Halfspace bot;        /**< 海底半空间参数 */
    Boundary2D top_curve; /**< 海面起伏曲线 (对应 .ati) */
    Boundary2D bot_curve; /**< 海底地形曲线 (对应 .bty) */
};

/** @struct SBP 对应 .sbp 文件内容 */
struct SBP {
    SBPFlag flag = SBPFlag::None; /**< 是否启用源方向图 */
    bool in_db   = false;         /**< angle_level_pairs 中的 level 是否为 dB 值 */
    std::vector<double> angle_level_pairs = {}; /**< {角度0, 幅度0, 角度1, 幅度1, ...} */
};

/** @struct ReflectionCoef 对应 .brc/.trc 文件中的单行数据 */
struct ReflectionCoef { double theta = 0.0; double r = 0.0; double phi = 0.0; };

/** @struct ReflectionTable 对应 .brc/.trc 文件内容 */
struct ReflectionTable { bool in_degrees = true; std::vector<ReflectionCoef> table = {}; };

/** @struct Reflection 对应 env 文件中引用 .brc/.trc 的设置 */
struct Reflection { ReflectionTable top; ReflectionTable bot; };

/** @struct Input2D 2D 计算的完整输入参数 */
struct Input2D {
    // ------------ env 兼容字段 -------------
    int n_media = 1;               /**< NMEDIA */
    /** @brief 对应 env 文件的 OPTIONS1 字符串, 如 "QVMT" */
    std::string options = "";

    /** @brief 若 angles.alpha 为空，则用 n_rays/start/end 自动生成等间隔射线角 */
    int n_rays = 0;
    double start_angle_deg = -90.0;
    double end_angle_deg   = 90.0;

    // ------------- 结构化主要字段 -------------
    Title title;
    Freq0 freq0;
    RunType run_type;
    Beam beam;
    Positions2D pos;
    Angles angles;
    SSP1D ssp;        // 即使是 Type 'Q'，也需要用这里的深度信息
    SSPQuad ssp_quad; // 仅当 options 中包含 'Q' 时使用
    Boundaries2D boundaries;
    SBP sbp;
    Reflection reflection;
};

// -----------------------------
// 输出
// -----------------------------

/** @struct TLResult2D TL 计算结果 */
struct TLResult2D {
    int width = 0;  /**< 距离方向的点数 */
    int height = 0; /**< 深度方向的点数 */
    std::vector<float> tl_db; /**< TL值 (dB)，大小为 width*height，布局: tl_db[iz*width + ir] */

    /** @brief 安全地访问指定点的 TL 值 */
    float at(int ir, int iz) const {
        return tl_db.at(static_cast<size_t>(iz) * static_cast<size_t>(width) + static_cast<size_t>(ir));
    }
};

// -----------------------------
// 计算接口
// -----------------------------

/**
 * @brief 执行 2D TL 计算
 * @param in 包含所有环境和计算参数的输入结构体
 * @return 包含 TL 声场结果的 TLResult2D 对象
 * @throw bhc_cpp::Error 如果计算失败，将抛出异常，异常消息包含 Bellhop 内部的详细错误信息
 */
BHC_CPP_API TLResult2D compute_tl_2d(const Input2D &in);

} // namespace bhc_cpp
