#include <cmath>
#include <iostream>
#include <vector>

#include "bhc/bhc_cpp_api.hpp"
#include "tl_png.hpp"

#ifdef _WIN32
#include <Windows.h>
#endif

// 深水 Q-SSP 示例（range 用 km，深度用 m）
// - Q-SSP: 完整、平滑的 c(z,r) 矩阵（避免外推/常数填充造成图像无结构）
// - 海底：非平直（起伏），并在 PNG 中以灰色遮罩显示

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    using namespace bhc_cpp;

    Input2D in;
    in.title.text = "Advanced(Q-SSP) Example";
    in.freq0.hz = 1000.0;

    // 启用 Q 型 SSP
    in.options = "Q";

    in.run_type.use_raw = false;
    in.run_type.e.tl = TLCoherence::Incoherent;
    in.run_type.e.infl = Influence::GeomGaussianCartesian;

    // box：range 用 km（10km），depth 用 m（5200m）
    in.beam.box_x_m = 10.0;
    in.beam.box_z_m = 5200.0;

    // 角度：自动生成
    in.n_rays = 121;
    in.start_angle_deg = -35.0;
    in.end_angle_deg = 35.0;

    // 源/接收（range: km，depth: m）
    in.pos.sz_m = {50.0f};

    in.pos.rr_in_km = true;
    in.pos.rr_m.resize(1001);
    for(int i = 0; i < 1001; ++i) in.pos.rr_m[i] = float(i) * 0.01f; // 0..10 km

    in.pos.rz_m.resize(501);
    for(int i = 0; i < 501; ++i) in.pos.rz_m[i] = float(i) * 10.0f; // 0..5000 m

    // -----------------------------
    // Q-SSP ranges（km）
    // -----------------------------
    in.ssp_quad.ranges_in_km = true;
    in.ssp_quad.ranges = {
        -115.61549999999983,
        0.0,
        1.9999999999999969,
        3.9999999999999938,
        5.9999999999999911,
        7.9999999999999876,
        9.999999999999984,
        11.999999999999982,
        13.999999999999979,
        15.999999999999975,
        17.999999999999972,
        19.999999999999968,
        21.999999999999964,
        115.61549999999983,
    };

    const int n_ranges = static_cast<int>(in.ssp_quad.ranges.size());

    // -----------------------------
    // 深度节点（m）
    // -----------------------------
    const std::vector<double> z_nodes = {
        0, 1, 2, 3, 5, 6, 7, 9, 11, 13, 15, 18, 21, 25, 29, 34, 40, 47, 55, 65, 77, 92,
        109, 130, 155, 186, 222, 266, 318, 380, 453, 541, 643, 763, 902, 1062, 1245, 1452,
        1684, 1941, 2225, 2533, 2865, 3220, 3597, 3992, 4405, 4833, 4891,
    };
    const int n_depths = static_cast<int>(z_nodes.size());

    // 生成 1D 基准剖面 c0(z)：深水典型 U 型
    auto c0_of_z = [](double z) {
        const double c_surface = 1546.0;
        const double c_min = 1482.0;
        const double z_min = 1000.0;
        const double scale = 3500.0;
        double dz = (z - z_min) / scale;
        return c_min + (c_surface - c_min) * (dz * dz + 0.15);
    };

    // range 扰动幅度随深度衰减
    auto amp_of_z = [](double z) {
        return 2.0 * std::exp(-z / 800.0) + 0.3; // m/s
    };

    constexpr double kPi = 3.141592653589793238462643383279502884;
    const double L_km = 21.0;

    in.ssp.type = SSPType::Quad;
    in.ssp.pts.resize(z_nodes.size());
    in.ssp_quad.ssp_matrix.resize(static_cast<size_t>(n_depths) * n_ranges);

    for(int iz = 0; iz < n_depths; ++iz) {
        double z = z_nodes[iz];
        double c0 = c0_of_z(z);
        double A  = amp_of_z(z);

        in.ssp.pts[static_cast<size_t>(iz)].z_m = z;
        in.ssp.pts[static_cast<size_t>(iz)].alphaR_mps = c0;

        for(int ir = 0; ir < n_ranges; ++ir) {
            double r_km = in.ssp_quad.ranges[static_cast<size_t>(ir)];
            double phase = 2.0 * kPi * (r_km / L_km);
            double c = c0 + A * std::sin(phase);
            in.ssp_quad.ssp_matrix[static_cast<size_t>(iz) * n_ranges + ir] = c;
        }
    }

    // -----------------------------
    // 边界：平直海面 + 非平直海底（起伏）
    // -----------------------------
    in.boundaries.top.bc = BoundaryCondition::Vacuum;
    in.boundaries.bot.bc = BoundaryCondition::Rigid;

    in.boundaries.top_curve.interp = BoundaryInterp::Linear;
    in.boundaries.top_curve.flag = BoundaryFlag::None;
    in.boundaries.top_curve.range_in_km = true;
    in.boundaries.top_curve.r = {0.0, 10.0};
    in.boundaries.top_curve.z = {0.0, 0.0};

    in.boundaries.bot_curve.interp = BoundaryInterp::Linear;
    in.boundaries.bot_curve.flag = BoundaryFlag::None;
    in.boundaries.bot_curve.range_in_km = true;

    // 50 个起伏点（0..10km），围绕 4800m 起伏
    const int Nb = 50;
    in.boundaries.bot_curve.r.resize(Nb);
    in.boundaries.bot_curve.z.resize(Nb);

    const double base_bot = 4800.0;
    for(int i = 0; i < Nb; ++i) {
        double t = (Nb == 1) ? 0.0 : double(i) / double(Nb - 1);
        double r_km = 10.0 * t;
        double und = 80.0 * std::sin(2.0 * kPi * t) + 40.0 * std::sin(2.0 * kPi * 3.0 * t + 0.7);
        double z = base_bot + und;
        in.boundaries.bot_curve.r[i] = r_km;
        in.boundaries.bot_curve.z[i] = z;
    }

    // 确保 SSP 最深点覆盖海底最大深度（给一点裕量）
    // 如果海底比 SSP 深，会在 preprocess 报 Boundary drops below SSP。
    double max_bot = in.boundaries.bot_curve.z[0];
    for(double z : in.boundaries.bot_curve.z) max_bot = max(max_bot, z);
    if(!in.ssp.pts.empty() && in.ssp.pts.back().z_m < max_bot + 50.0) {
        // 追加一个深度点，声速沿用最后一层的 c0
        SSP1DPoint p = in.ssp.pts.back();
        p.z_m = max_bot + 50.0;
        in.ssp.pts.push_back(p);

        // 同步扩展 Q-SSP 矩阵：为新增深度行填同样的 c0(z)
        const int new_depths = static_cast<int>(in.ssp.pts.size());
        std::vector<double> new_mat(static_cast<size_t>(new_depths) * n_ranges);
        // 复制旧数据
        for(int iz = 0; iz < n_depths; ++iz)
            for(int ir = 0; ir < n_ranges; ++ir)
                new_mat[static_cast<size_t>(iz) * n_ranges + ir] = in.ssp_quad.ssp_matrix[static_cast<size_t>(iz) * n_ranges + ir];
        // 新行
        double znew = p.z_m;
        double c0 = c0_of_z(znew);
        double A  = amp_of_z(znew);
        for(int ir = 0; ir < n_ranges; ++ir) {
            double rkm = in.ssp_quad.ranges[static_cast<size_t>(ir)];
            double phase = 2.0 * kPi * (rkm / L_km);
            new_mat[static_cast<size_t>(new_depths - 1) * n_ranges + ir] = c0 + A * std::sin(phase);
        }
        in.ssp_quad.ssp_matrix.swap(new_mat);
    }

    // -----------------------------
    // 运行 + 输出 PNG（含海底遮罩）
    // -----------------------------
    try {
        auto out = compute_tl_2d(in);
        std::cout << "Advanced Q-SSP example: calculation successful! Grid: "
                  << out.width << "x" << out.height << "\n";

        const std::string png_name = "example_advanced.png";
        if(!tl_png::write_tl_png(png_name, out.width, out.height, out.tl_db)) {
            std::cerr << "Failed to write PNG: " << png_name << "\n";
            return 2;
        }
        std::cout << "Wrote PNG: " << png_name << "\n";

    } catch(const std::exception &e) {
        std::cerr << "Advanced example failed: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
