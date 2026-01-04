#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "TEDSbellhop/TEDSbellhop.hpp"

#include "tl_png.hpp"

#ifdef _WIN32
#include <windows.h>
#endif

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    using namespace tedsbellhop;

    TL2DInput in;

    in.title = "TEDSbellhop TL 示例 (不读文件)";
    in.freq_hz = 2000.0;

    // 不需要 env 文件；options 按需传（这里不使用 Q）
    in.options = "";

    in.run_type.use_raw = false;
    in.run_type.tl = TLCoherence::Incoherent;
    in.run_type.infl = Influence::GeomGaussianCartesian;

    in.beam.box_x_m = 10000.0;
    in.beam.box_z_m = 500.0;

    in.pos.source_depths_m = {50.0f};

    in.pos.receiver_ranges.resize(501);
    for(int i = 0; i < 501; ++i) in.pos.receiver_ranges[i] = float(i) * 20.0f; // 0-10km, 20m
    in.pos.receiver_depths_m.resize(301);
    for(int i = 0; i < 301; ++i) in.pos.receiver_depths_m[i] = float(i) * 1.0f; // 0-300m

    //in.angles.alpha.resize(121);
    //for(int i = 0; i < 121; ++i) in.angles.alpha[i] = -35.0 + 70.0 * double(i) / 120.0;
    in.n_rays = 2000;
    in.start_angle_deg = -90.f;
    in.end_angle_deg = 90.f;

    in.ssp.type = SSPType::CLinear;
    in.ssp.points = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {50.0, 1510.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {200.0, 1490.0, 0.0, 1.0, 0.0, 0.0},
        {250.0, 1470.0, 0.0, 1.0, 0.0, 0.0},
        {300.0, 1490.0, 0.0, 1.0, 0.0, 0.0},
        {340.0, 1520.0, 0.0, 1.0, 0.0, 0.0},
        {400.0, 1530.0, 0.0, 1.0, 0.0, 0.0},
    };

    in.boundaries.top.bc = BoundaryCondition::Vacuum;
    in.boundaries.bottom.bc = BoundaryCondition::Rigid;

    in.boundaries.top_curve.interp = BoundaryInterp::Curvilinear;
    in.boundaries.top_curve.r = {0.0, 10000.0};
    in.boundaries.top_curve.z = {0.0, 0.0};

    //in.boundaries.bottom_curve.interp = BoundaryInterp::Linear;
    //in.boundaries.bottom_curve.r = {0.0, 4000.0, 10000.0};
    //in.boundaries.bottom_curve.z = {280.0, 300.0, 250.0};


    in.boundaries.bottom_curve.interp = BoundaryInterp::Linear;

    // 50 个起伏点（0..10000m），围绕 280m 上下起伏
    const int Nb = 500;
    in.boundaries.bottom_curve.r.resize(Nb);
    in.boundaries.bottom_curve.z.resize(Nb);

    const double base_bot = 280.0;                  // 平均海底深度
    const double A1 = 55.0;                         // 主起伏幅度
    const double A2 = 12.0;                         // 次起伏幅度
    constexpr double kPi = 3.14159265358979323846;

    for (int i = 0; i < Nb; ++i) {
        double t = (Nb == 1) ? 0.0 : double(i) / double(Nb - 1);
        double rr = 10000.0 * t;                    // 0..10000 m

        // 连续起伏：两种频率叠加 + 相位偏移（更“自然”）
        double und = A1 * std::sin(2.0 * kPi * t)
            + A2 * std::sin(2.0 * kPi * 3.0 * t + 0.7);

        in.boundaries.bottom_curve.r[i] = rr;
        in.boundaries.bottom_curve.z[i] = base_bot + und;
    }

    // 捕获底层输出（可选）
    in.prt_callback = [](const char* m) {
        if(m) std::cout << m;
    };
    in.output_callback = [](const char* m) {
        if(m) std::cout << "Out: " << m << "\n";
    };

    try {
        auto out = compute_tl_2d(in);
        std::cout << "计算成功! 网格尺寸: " << out.width << "x" << out.height << "\n";
        std::cout << "TL(5km,50m)=" << out.at(250, 50) << " dB\n";

        // 简单 sanity check
        if(out.tl_db.empty()) return 2;
        if(!std::isfinite(out.tl_db[out.tl_db.size()/2])) return 3;

        // 保存为图片（不叠加地形遮罩，直接按 TL 输出）
        const std::string png_name = "example_teds_tl_cpp.png";
        if(!tl_png::write_tl_png(png_name, out.width, out.height, out.tl_db)) {
            std::cerr << "写 PNG 失败: " << png_name << "\n";
            return 4;
        }
        std::cout << "已生成图片: " << png_name << "\n";

    } catch(const std::exception &e) {
        std::cerr << "计算失败: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

