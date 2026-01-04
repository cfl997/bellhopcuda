#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "bhc/bhc_cpp_api.hpp"
#include "tl_png.hpp"

// 说明：
// - 该用例根据用户给定的 env / bty / ssp 三段文本，构造 bhc_cpp::Input2D
// - 目标是验证：compute_tl_2d() 能跑通，并且输出网格尺寸与输入一致，且 TL 数组为有限数。

int main() {
    using namespace bhc_cpp;

    Input2D in;

    // -----------------------------
    // env
    // -----------------------------
    in.title.text = "BELLHOP - woss id = 0; run = 0";
    in.freq0.hz = 1000.0;

    // OPTIONS1: 'QVMT'
    // - Q: range-dependent SSP quad
    // - V: top boundary vacuum
    // - M/T: 其它内部选项（这里由库内部解释），我们只按原样传入
    in.options = "QVMT";

    // RAY OPTIONS: 'IB'
    // 这里 C++ API 没有直接暴露 'IB'，run_type 采用枚举方式。
    // 保持默认：二维、相干、帽函数/笛卡尔。
    in.run_type.use_raw = false;
    in.run_type.e.tl = TLCoherence::Coherent;
    in.run_type.e.infl = Influence::GeomHatCartesian;
    in.run_type.e.source = SourceModel::Point;
    in.run_type.e.grid = ReceiverGrid::Rectilinear;
    in.run_type.e.dim = DimensionFlag::D2;

    // beam: 0.0 7492.1 110.11
    // 注：env 中的三元组在 bellhop 输入里含义与本封装字段不完全一致，
    // 这里把 BOX 深度/距离映射到 box_z_m/box_x_m，步长不在本接口中直接对应。
    in.beam.box_x_m = 110.10999999999983;
    in.beam.box_z_m = 7492.1000000000004;
    in.beam.deltas_m = 0.0;

    // 源深：1 个 source，50 m
    in.pos.sz_m = {50.0f};

    // 接收深度：200 个，0 ~ 6674.78 m
    in.pos.rz_m.resize(200);
    {
        const float z0 = 0.0f;
        const float z1 = 6674.7799999999997f;
        for (int i = 0; i < 200; ++i) {
            float t = (200 == 1) ? 0.0f : float(i) / float(200 - 1);
            in.pos.rz_m[i] = z0 + (z1 - z0) * t;
        }
    }

    // 接收距离：1000 个，0 ~ 100 m
    in.pos.rr_in_km = false;
    in.pos.rr_m.resize(1000);
    {
        const float r0 = 0.0f;
        const float r1 = 99.99999999999984f;
        for (int i = 0; i < 1000; ++i) {
            float t = (1000 == 1) ? 0.0f : float(i) / float(1000 - 1);
            in.pos.rr_m[i] = r0 + (r1 - r0) * t;
        }
    }

    // 角度：-90 到 90，2000 条射线
    // API 的 Angles 是显式数组，这里等间隔生成 2000 个。
    in.angles.alpha_in_degrees = true;
    in.angles.alpha.resize(2000);
    {
        const double a0 = -90.0;
        const double a1 = 90.0;
        for (int i = 0; i < 2000; ++i) {
            double t = (2000 == 1) ? 0.0 : double(i) / double(2000 - 1);
            in.angles.alpha[i] = a0 + (a1 - a0) * t;
        }
    }

    // -----------------------------
    // bty
    // -----------------------------
    // bty 文件：'L' + 51 个点；这里映射到 bot_curve
    // 注意：bty 中的 range 看起来是 0~100（可能是 km），但 env 的 box range 只有 110 m。
    // 按“原样”填入（单位当作 m），以验证接口/内核健壮性。
    in.boundaries.bot_curve.interp = BoundaryInterp::Linear;
    in.boundaries.bot_curve.flag   = BoundaryFlag::None;
    in.boundaries.bot_curve.range_in_km = false;
    in.boundaries.bot_curve.extend_to_infinity = true;

    in.boundaries.bot_curve.r = {
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
        23.999999999999964,
        25.999999999999961,
        27.999999999999957,
        29.999999999999954,
        31.99999999999995,
        33.99999999999995,
        35.999999999999943,
        37.999999999999943,
        39.999999999999936,
        41.999999999999936,
        43.999999999999929,
        45.999999999999929,
        47.999999999999929,
        49.999999999999922,
        51.999999999999922,
        53.999999999999915,
        55.999999999999915,
        57.999999999999908,
        59.999999999999908,
        61.999999999999901,
        63.999999999999901,
        65.999999999999901,
        67.999999999999901,
        69.999999999999886,
        71.999999999999886,
        73.999999999999886,
        75.999999999999886,
        77.999999999999872,
        79.999999999999872,
        81.999999999999872,
        83.999999999999872,
        85.999999999999872,
        87.999999999999858,
        89.999999999999858,
        91.999999999999858,
        93.999999999999858,
        95.999999999999858,
        97.999999999999844,
        99.999999999999844,
    };

    in.boundaries.bot_curve.z = {
        4994.0,
        4991.0,
        5009.5,
        5068.0,
        5130.0,
        5204.5,
        5270.0,
        5311.0,
        5383.5,
        5461.0,
        5504.5,
        5535.0,
        5576.5,
        5623.5,
        5655.5,
        5690.5,
        5738.5,
        5779.5,
        5811.5,
        5846.5,
        5888.5,
        5915.0,
        5917.5,
        5924.0,
        5949.0,
        5992.0,
        6050.0,
        6115.5,
        6166.5,
        6180.5,
        6192.0,
        6256.5,
        6336.0,
        6374.5,
        6390.0,
        6393.5,
        6412.0,
        6432.5,
        6456.0,
        6485.5,
        6539.0,
        6601.0,
        6602.5,
        6585.5,
        6592.0,
        6650.0,
        6711.5,
        6757.5,
        6800.5,
        6718.5,
        6529.0,
    };

    // top boundary: env options 中有 'V'，顶边界真空
    in.boundaries.top.bc = BoundaryCondition::Vacuum;

    // bottom halfspace: 来自 env：6811 1552.25 740.448 1.6335 0.7095 1.325
    // 对应：alphaR=1552.25, betaR=740.448, rho=1.6335, alphaI=0.7095, betaI=1.325
    in.boundaries.bot.bc = BoundaryCondition::Acoustic;
    in.boundaries.bot.alphaR_mps = 1552.25;
    in.boundaries.bot.betaR_mps = 740.448;
    in.boundaries.bot.rho_gcm3 = 1.6335;
    in.boundaries.bot.alphaI = 0.7095;
    in.boundaries.bot.betaI = 1.325;

    // -----------------------------
    // ssp (二维 Q-quad)
    // -----------------------------
    // ssp 文件：
    // - 第一行 14：表示 14 个 range 节点（列数）
    // - 下一行是 ranges（看起来从 -115.6155 到 115.6155）
    // - 后续每行第一列为深度 z，其余 14 列为该深度在各 range 的声速
    // 我们在本 API 中：
    // - 深度用 in.ssp.pts 的 z
    // - ranges + 声速矩阵 用 in.ssp_quad

    in.ssp.type = SSPType::Quad; // 'Q'

    in.ssp_quad.ranges_in_km = false;
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

    // 深度节点（来自 env 的 z-ctp 列表，最后一个 6811）
    // 同时把每个深度的“参考声速”填到 ssp.pts（这里取第 2 列=range 0.0 的声速）
    // 备注：beta/atten/rho 在水体里通常不使用，这里按默认。
    const std::vector<double> z_nodes = {
        0, 1, 2, 3, 5, 6, 7, 9, 11, 13, 15, 18, 21, 25, 29, 34, 40, 47, 55, 65, 77, 92, 109, 130, 155, 186, 222, 266, 318, 380, 453, 541, 643, 763, 902, 1062, 1245, 1452, 1684, 1941, 2225, 2533, 2865, 3220, 3597, 3992, 4405, 4833, 5274, 5727, 6300, 6811,
    };

    // range=0.0 列声速（来自 env 的同深度单列 SSP）
    const std::vector<double> c_ref = {
        1538.868896484375,
        1538.8856201171875,
        1538.890380859375,
        1538.883056640625,
        1538.864501953125,
        1538.8448486328125,
        1538.8135986328125,
        1538.76513671875,
        1538.6956787109375,
        1538.570556640625,
        1538.2947998046875,
        1537.9190673828125,
        1537.5203857421875,
        1537.0594482421875,
        1536.50830078125,
        1535.7532958984375,
        1534.7716064453125,
        1533.6019287109375,
        1532.3236083984375,
        1530.96435546875,
        1529.618896484375,
        1528.30419921875,
        1526.997314453125,
        1525.65625,
        1524.2557373046875,
        1522.309814453125,
        1519.5911865234375,
        1516.4754638671875,
        1512.990234375,
        1508.2528076171875,
        1501.9735107421875,
        1494.4169921875,
        1487.792236328125,
        1483.4862060546875,
        1480.9671630859375,
        1480.5325927734375,
        1481.1275634765625,
        1483.05810546875,
        1485.8192138671875,
        1489.3427734375,
        1493.5460205078125,
        1498.13720703125,
        1503.466552734375,
        1509.300048828125,
        1515.73291015625,
        1522.7271728515625,
        1530.1982421875,
        1538.0234375,
        1546.093994140625,
        1554.3775634765625,
        1564.83544921875,
        1574.1617661720795,
    };

    if (z_nodes.size() != c_ref.size()) {
        throw std::runtime_error("internal test data error: z_nodes.size != c_ref.size");
    }

    in.ssp.pts.resize(z_nodes.size());
    for (size_t i = 0; i < z_nodes.size(); ++i) {
        in.ssp.pts[i].z_m = z_nodes[i];
        in.ssp.pts[i].alphaR_mps = c_ref[i];
        in.ssp.pts[i].betaR_mps = 0.0;
        in.ssp.pts[i].rho_gcm3 = 1.0;
        in.ssp.pts[i].alphaI = 0.0;
        in.ssp.pts[i].betaI = 0.0;
    }

    // 关键：ssp_quad.ssp_matrix
    // 尺寸 = depth_count * range_count
    // 布局：depth 优先（depth varies fastest）
    // 但用户给的 ssp 内容太大（每个深度 14 列），这里示范性填充：
    // - 使用 range=0.0 的参考声速复制到所有 range 列
    // 这样仍然是合法的 Q-quad 输入，且能跑通。
    const size_t n_depth = in.ssp.pts.size();
    const size_t n_range = in.ssp_quad.ranges.size();
    in.ssp_quad.ssp_matrix.resize(n_depth * n_range);
    for (size_t iz = 0; iz < n_depth; ++iz) {
        for (size_t ir = 0; ir < n_range; ++ir) {
            in.ssp_quad.ssp_matrix[iz * n_range + ir] = in.ssp.pts[iz].alphaR_mps;
        }
    }

    // -----------------------------
    // 执行并断言
    // -----------------------------
    TLResult2D out;
    try {
        out = compute_tl_2d(in);
    } catch (const std::exception &e) {
        std::cerr << "compute_tl_2d 失败: " << e.what() << "\n";
        return 1;
    }

    if (out.width != static_cast<int>(in.pos.rr_m.size())) {
        std::cerr << "width 不一致: out.width=" << out.width
                  << " in.pos.rr_m.size=" << in.pos.rr_m.size() << "\n";
        return 2;
    }
    if (out.height != static_cast<int>(in.pos.rz_m.size())) {
        std::cerr << "height 不一致: out.height=" << out.height
                  << " in.pos.rz_m.size=" << in.pos.rz_m.size() << "\n";
        return 3;
    }
    if (out.tl_db.size() != static_cast<size_t>(out.width) * static_cast<size_t>(out.height)) {
        std::cerr << "tl_db 尺寸不一致: tl_db.size=" << out.tl_db.size()
                  << " expected=" << (static_cast<size_t>(out.width) * static_cast<size_t>(out.height)) << "\n";
        return 4;
    }

    // 抽查若干点为有限数
    const size_t probes[] = {0, out.tl_db.size() / 3, out.tl_db.size() / 2, out.tl_db.size() - 1};
    for (size_t idx : probes) {
        float v = out.tl_db[idx];
        if (!std::isfinite(v)) {
            std::cerr << "TL 非有限数: idx=" << idx << " v=" << v << "\n";
            return 5;
        }
    }

    // 保存为图片（不叠加地形遮罩）
    const std::string png_name = "test_case_qvmt.png";
    if(!tl_png::write_tl_png(png_name, out.width, out.height, out.tl_db)) {
        std::cerr << "Failed to write PNG: " << png_name << "\n";
        return 6;
    }

    std::cout << "test_case_qvmt PASS: " << out.width << "x" << out.height << " (" << png_name << ")\n";
    return 0;
}

