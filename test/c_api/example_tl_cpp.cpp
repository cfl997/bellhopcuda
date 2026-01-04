#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "bhc/bhc_cpp_api.hpp"

#ifdef _WIN32
#include <windows.h>
#endif

#include "tl_png.hpp"

// Windows 宏污染修复：windows.h 可能定义 min/max，导致 std::min/std::max 解析异常
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif


int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    using namespace bhc_cpp;

    Input2D in;

    in.title.text = "C++ API TL 示例 (含海底地形)";
    in.freq0.hz = 1000.0;

    in.run_type.use_raw = false;
    in.run_type.e.tl = TLCoherence::Incoherent;
    in.run_type.e.infl = Influence::GeomGaussianCartesian;

    in.beam.box_x_m = 10000.0;
    in.beam.box_z_m = 320.0; // 略大于最大深度

    in.pos.sz_m = {50.0f};
    in.pos.rr_m.resize(501);
    for(int i=0;i<501;++i) in.pos.rr_m[i] = float(i) * 20.0f; // 0-10km, 20m步长
    in.pos.rz_m.resize(301);
    for(int i=0;i<301;++i) in.pos.rz_m[i] = float(i) * 1.0f; // 0-300m, 1m步长

    in.angles.alpha.resize(121);
    for(int i=0;i<121;++i) in.angles.alpha[i] = -35.0 + 70.0 * double(i) / 120.0;

    in.ssp.type = SSPType::CLinear;
    in.ssp.pts = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1495.0, 0.0, 1.0, 0.0, 0.0},
        {250.0, 1485.0, 0.0, 1.0, 0.0, 0.0},
        {300.0, 1480.0, 0.0, 1.0, 0.0, 0.0},
    };

    in.boundaries.top.bc = BoundaryCondition::Vacuum;
    in.boundaries.bot.bc = BoundaryCondition::Rigid;

    in.boundaries.top_curve.interp = BoundaryInterp::Linear;
    in.boundaries.top_curve.flag   = BoundaryFlag::None;
    in.boundaries.top_curve.r = {0.0, 10000.0};
    in.boundaries.top_curve.z = {0.0, 0.0};

    // 新增：一个简单的上坡地形
    in.boundaries.bot_curve.interp = BoundaryInterp::Linear;
    in.boundaries.bot_curve.flag   = BoundaryFlag::None;
    in.boundaries.bot_curve.r = {0.0, 4000.0, 10000.0};
    in.boundaries.bot_curve.z = {280.0, 300.0, 250.0};

    try {
        auto out = compute_tl_2d(in);
        std::cout << "计算成功! 网格尺寸: " << out.width << "x" << out.height << "\n";
        std::cout << "TL(5km,50m)=" << out.at(250, 50) << " dB\n";

        const std::string img_path = "example_tl_cpp.png";
        if(!tl_png::write_tl_png(img_path, out.width, out.height, out.tl_db)) {
            std::cerr << "写 PNG 失败: " << img_path << "\n";
            return 1;
        }


        std::cout << "已生成图片: " << img_path << "\n";

    } catch(const std::exception &e) {
        std::cerr << "计算失败: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
