#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <string>
#include <filesystem>

#include "TEDSbellhop/TEDSbellhop.hpp"
#include "tl_png.hpp"

#ifdef _WIN32
#include <windows.h>
#endif

// 辅助函数：创建一个基础的输入配置
tedsbellhop::TL2DInput create_base_input() {
    using namespace tedsbellhop;
    TL2DInput in;

    in.title = "TEDSbellhop Parallel Test";
    in.freq_hz = 1000.0;
    in.options = "";

    in.run_type.use_raw = false;
    in.run_type.tl = TLCoherence::Incoherent;
    in.run_type.infl = Influence::GeomGaussianCartesian;

    in.beam.box_x_m = 10000.0;
    in.beam.box_z_m = 500.0;

    in.pos.receiver_ranges.resize(501);
    for (int i = 0; i < 501; ++i) in.pos.receiver_ranges[i] = float(i) * 20.0f; // 0-10km, 20m
    in.pos.receiver_depths_m.resize(301);
    for (int i = 0; i < 301; ++i) in.pos.receiver_depths_m[i] = float(i) * 1.0f; // 0-300m

    //in.angles.alpha.resize(121);
    //for(int i = 0; i < 121; ++i) in.angles.alpha[i] = -35.0 + 70.0 * double(i) / 120.0;

    in.n_rays = 2000;
    in.start_angle_deg = -90.f;
    in.end_angle_deg = 90.f;

    in.ssp.type = SSPType::CLinear;
    in.ssp.points = {
        {0.0,   1500.0, 0.0, 1.0, 0.0, 0.0},
        {50.0,  1510.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {200.0, 1490.0, 0.0, 1.0, 0.0, 0.0},
        {250.0, 1470.0, 0.0, 1.0, 0.0, 0.0},
        {300.0, 1490.0, 0.0, 1.0, 0.0, 0.0},
        {340.0, 1520.0, 0.0, 1.0, 0.0, 0.0},
        {400.0, 1530.0, 0.0, 1.0, 0.0, 0.0},
    };

    in.boundaries.top.bc = BoundaryCondition::Vacuum;
    in.boundaries.bottom.bc = BoundaryCondition::Rigid;

    in.boundaries.top_curve.r = {0.0, 10000.0};
    in.boundaries.top_curve.z = {0.0, 0.0};

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

    return in;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    const int num_profiles = 8;

    // 输出目录：test_teds_parallel_outputs
    const std::filesystem::path out_dir = "test_teds_parallel_outputs";
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if(ec) {
        std::cerr << "创建输出目录失败: " << out_dir.string() << " : " << ec.message() << "\n";
        return 2;
    }

    std::vector<tedsbellhop::TL2DJob> jobs;
    jobs.reserve(num_profiles);

    std::cout << "启动 " << num_profiles << " 个并行计算任务...\n";

    for (int i = 0; i < num_profiles; ++i) {
        auto in = create_base_input();

        // 为每个任务创建一个独特的输入（示例：改变声源深度）
        in.pos.source_depths_m = {static_cast<float>(20 + i * 5)};
        in.title = "TEDSbellhop Parallel Test - Profile " + std::to_string(i);

        std::cout << "  - 启动剖面 " << i << " (声源深度: " << in.pos.source_depths_m[0] << "m)\n";
        jobs.push_back(tedsbellhop::start_tl_2d(in));
    }

    std::cout << "\n所有任务已启动，开始轮询进度：\n";

    int completed_count = 0;
    std::vector<bool> job_done(num_profiles, false);

    while (completed_count < num_profiles) {
        std::string progress_line = "进度: ";
        completed_count = 0;

        for (int i = 0; i < num_profiles; ++i) {
            if (!job_done[i]) {
                if (jobs[i].ready()) {
                    job_done[i] = true;
                }
            }
            if(job_done[i]) completed_count++;

            int p = jobs[i].progress();
            progress_line += "P" + std::to_string(i) + ":";
            if (p < 0) progress_line += "--";
            else if (p < 10) progress_line += "  " + std::to_string(p);
            else if (p < 100) progress_line += " " + std::to_string(p);
            else progress_line += std::to_string(p);
            progress_line += "% | ";
        }

        std::cout << "\r" << progress_line << std::flush;
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    std::cout << "\n\n所有任务已完成。正在获取结果并保存图片...\n";

    bool all_ok = true;
    for (int i = 0; i < num_profiles; ++i) {
        try {
            auto result = jobs[i].get();
            std::cout << "  - 剖面 " << i << " 成功: 网格 " << result.width << "x" << result.height << "\n";

            const std::filesystem::path png_path = out_dir / ("profile_" + std::to_string(i) + ".png");
            auto se = result.tl_db;
            //for (auto& data : se)
            //    data = 80 - data;
            if(!tl_png::write_tl_png(png_path.string(), result.width, result.height, se)) {
                std::cerr << "    保存 PNG 失败: " << png_path.string() << "\n";
                all_ok = false;
            } else {
                std::cout << "    已保存: " << png_path.string() << "\n";
            }

        } catch (const std::exception& e) {
            std::cerr << "  - 剖面 " << i << " 失败: " << e.what() << "\n";
            all_ok = false;
        }
    }

    if (all_ok) {
        std::cout << "\n所有剖面均计算成功！输出目录: " << out_dir.string() << "\n";
        return 0;
    } else {
        std::cerr << "\n部分剖面计算/保存失败。输出目录: " << out_dir.string() << "\n";
        return 1;
    }
}
