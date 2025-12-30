#pragma once

// 说明：这是 bellhopcxx/bellhopcuda 的 C++17 友好封装接口。
// 目标：完全不读任何输入文件（.env/.bty/.ati/.ssp/.brc/.trc/.sbp...），
//      外部以结构体+vector/string 形式提供全部参数，库内部拷贝并写入 bhcParams。

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------
// DLL 导出/导入宏（Windows 必需，否则 exe 链接会 LNK2019）
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
// 基础工具
// -----------------------------

struct Error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

// -----------------------------
// 与 "文件语义" 等价的输入结构（全部由调用者显式填写，struct 里给默认值）
// -----------------------------

// env: Title
struct Title {
    std::string text = ""; // 将截断/填充到 bellhop 的 80 字符 Title
};

// env: Freq0
struct Freq0 {
    double hz = 1000.0;
};

// env: RunType（7字符，保持与原 bellhop env 完全一致）
struct RunType {
    // 库内部会取前 7 字符，不足补空格
    std::string value = "C"; // 默认：相干 TL
};

// env: Beam/Box 等
struct Beam {
    // Beam->Type[4]，保持与原 env 语义一致，不足补空格
    std::string type = "G";

    // 射线/场计算范围框（单位：米）
    double box_x_m = 10000.0;
    double box_z_m = 200.0;

    // 步长（单位：米），<=0 表示自动选择
    double deltas_m = -1.0;

    // 其它高级参数（按需扩展）
    double eps_multiplier = 1.0;
};

// env: 源/接收器位置（2D）
struct Positions2D {
    std::vector<float> sz_m = {50.0f};
    std::vector<float> rr_m = {}; // 必填
    std::vector<float> rz_m = {}; // 必填

    // 对应 env 的 RrInKm
    bool rr_in_km = false;
};

// env: RayAngles（这里只先实现 elevation alpha）
struct Angles {
    std::vector<double> alpha = {}; // 必填
    bool alpha_in_degrees = true;
};

// env: SSP 1D 点（与 env 中 SSP 行语义一致）
struct SSP1DPoint {
    double z_m = 0.0;
    double alphaR_mps = 1500.0;
    double betaR_mps  = 0.0;
    double rho_gcm3   = 1.0;
    double alphaI     = 0.0;
    double betaI      = 0.0;
};

// env: SSP 1D（对应 env 中的主 SSP 段）
struct SSP1D {
    // ssp->Type: 'N','C','S','P','A' 等
    char type = 'C';

    // ssp->AttenUnit[2]：例如 "W " 表示 dB/wavelength
    std::string atten_unit = "W ";

    // 必填：至少 2 点，z 单调递增
    std::vector<SSP1DPoint> pts = {};
};

// bty/ati: 2D 边界曲线（与 .bty/.ati 文件语义一致）
// 规则（2D）：
// - type 必须长度>=1；type[0] 必须是 'C' 或 'L'（'R' 在 2D 不支持）；
// - type[1] 必须是 'S' / ' ' / 'L'；
// - r 必须单调递增；r/z 长度一致且 >=2；
// - 若 type[1]=='L'，则每个点必须提供半空间参数（hs.size()==N）。
struct Boundary2D {
    std::string type = "LS";
    bool range_in_km = false;

    // 点列（不包含两端扩展点）
    std::vector<double> r = {};
    std::vector<double> z = {};

    // 仅当 type[1]=='L' 时使用：每个点的半空间参数
    std::vector<SSP1DPoint> hs = {};

    // 是否按 bellhop 文件语义自动扩展到“±无穷”（默认 true）
    bool extend_to_infinity = true;

    // 扩展范围（当 extend_to_infinity=true 时生效）
    // 注意：如果 range_in_km=false，则单位为 m；如果 range_in_km=true，则单位为 km。
    double extend_left = 1e9;  // 约等于 -∞
    double extend_right = 1e9; // 约等于 +∞
};

// 顶/底边界条件（对应 env 的 BoundaryCondTop/Bottom + BotOpt/TopOpt 组合）
struct Halfspace {
    char bc = 'R';

    // bc=='A' 时使用：
    double alphaR_mps = 1600.0;
    double betaR_mps  = 0.0;
    double rho_gcm3   = 1.8;
    double alphaI     = 0.0;
    double betaI      = 0.0;
};

struct Boundaries2D {
    Halfspace top;
    Halfspace bot;

    // 额外曲线：top_curve 对应 .ati，bot_curve 对应 .bty
    Boundary2D top_curve;
    Boundary2D bot_curve;
};

// SBP: 源指向性（对应 .sbp 内容）
struct SBP {
    char flag = 'N';
    bool in_db = false;
    std::vector<double> angle_level_pairs = {}; // {angle0, level0, angle1, level1, ...}
};

// TRC/BRC 反射系数（对应 .trc/.brc 文件）
struct ReflectionCoef {
    double theta = 0.0;
    double r = 0.0;
    double phi = 0.0;
};

struct ReflectionTable {
    bool in_degrees = true;
    std::vector<ReflectionCoef> table = {};
};

struct Reflection {
    ReflectionTable top;
    ReflectionTable bot;
};

// 总输入（2D TL）
struct Input2D {
    Title title;
    Freq0 freq0;
    RunType run_type;
    Beam beam;
    Positions2D pos;
    Angles angles;
    SSP1D ssp;
    Boundaries2D boundaries;
    SBP sbp;
    Reflection reflection;
};

// -----------------------------
// 输出
// -----------------------------

struct TLResult2D {
    int width = 0;
    int height = 0;
    std::vector<float> tl_db;

    float at(int ir, int iz) const {
        return tl_db.at(static_cast<size_t>(iz) * static_cast<size_t>(width)
                        + static_cast<size_t>(ir));
    }
};

// -----------------------------
// 计算接口
// -----------------------------

BHC_CPP_API TLResult2D compute_tl_2d(const Input2D &in);

} // namespace bhc_cpp
