#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// stb_image_write：仅用于把 TL 结果写成 PNG
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Windows 宏污染修复：windows.h 可能定义 min/max，导致 std::min/std::max 解析异常
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

namespace tl_png {

struct RGB {
    uint8_t r, g, b;
};

static inline float clamp01(float x) { return std::max(0.0f, std::min(1.0f, x)); }

static inline RGB lerp(const RGB &a, const RGB &b, float t) {
    t = clamp01(t);
    return RGB{
        static_cast<uint8_t>(a.r + (b.r - a.r) * t),
        static_cast<uint8_t>(a.g + (b.g - a.g) * t),
        static_cast<uint8_t>(a.b + (b.b - a.b) * t),
    };
}

// =============================================
// 色卡：把下面任意一个 stops_* 复制到 colormap() 里就能切换
// 最直观手动方式：
// 1) 在 colormap() 中把 `const auto &stops = stops_A;` 改成 stops_B/C/... 即可。
// 2) 如需更多色卡，照这个格式再加一个 stops_X。
// =============================================

// A：你当前的 stops（偏“红-黄-绿-青-蓝”）
static const std::vector<std::pair<float, RGB>> stops_A = {
    {0.00f, {128,   0,   0}},
    {0.20f, {255,   0,   0}},
    {0.35f, {255, 215, 107}},
    {0.39f, {  0, 255,   0}},
    {0.48f, {  0, 230, 255}},
    {0.71f, {  0,   0, 255}},
    {1.00f, { 10,  20, 240}},
};

// B：高对比（黑蓝→紫→洋红→橙→黄→白）
static const std::vector<std::pair<float, RGB>> stops_B = {
    {0.00f, {  0,   0,   4}},
    {0.15f, { 40,  10, 107}},
    {0.40f, {148,  33, 107}},
    {0.65f, {237,  87,  59}},
    {0.85f, {253, 201, 100}},
    {1.00f, {255, 255, 255}},
};

// C：灰度（黑→白）
static const std::vector<std::pair<float, RGB>> stops_C = {
    {0.00f, {  0,   0,   0}},
    {1.00f, {255, 255, 255}},
};

// D：turbo-like（感知更均匀一些）
static const std::vector<std::pair<float, RGB>> stops_D = {
    {0.00f, { 48,  18,  59}},
    {0.20f, { 50,  70, 206}},
    {0.40f, { 36, 188, 235}},
    {0.60f, {149, 249,  74}},
    {0.80f, {253, 174,  35}},
    {1.00f, {180,   4,  38}},
};

// E：冷色系（黑→深蓝→青→白）适合强调弱信号/细节
static const std::vector<std::pair<float, RGB>> stops_E = {
    {0.00f, {  0,   0,   0}},
    {0.30f, {  0,   0, 120}},
    {0.55f, {  0, 180, 255}},
    {0.75f, {150, 240, 255}},
    {1.00f, {255, 255, 255}},
};

// F：暖色系（黑→紫→红→橙→黄→白）适合强调强信号/热点
static const std::vector<std::pair<float, RGB>> stops_F = {
    {0.00f, {  0,   0,   0}},
    {0.20f, { 60,   0, 100}},
    {0.45f, {200,   0,  60}},
    {0.70f, {255, 120,   0}},
    {0.85f, {255, 220,  80}},
    {1.00f, {255, 255, 255}},
};

// G：对称红蓝（中性点清晰，适合看“正负差值/残差”类图）
static const std::vector<std::pair<float, RGB>> stops_G = {
    {0.00f, {  0,   0, 180}},
    {0.50f, {240, 240, 240}},
    {1.00f, {180,   0,   0}},
};

static inline RGB colormap_from_stops(const std::vector<std::pair<float, RGB>> &stops, float t) {
    t = clamp01(t);
    for(size_t i = 1; i < stops.size(); ++i) {
        if(t <= stops[i].first) {
            float t0 = stops[i - 1].first;
            float t1 = stops[i].first;
            float u  = (t - t0) / (t1 - t0);
            return lerp(stops[i - 1].second, stops[i].second, u);
        }
    }
    return stops.back().second;
}

static inline RGB colormap(float t) {
    // ================================
    // 手动切换色卡只需要改这一行：
    // const auto &stops = stops_A;
    // 可选：stops_A / stops_B / stops_C / stops_D / stops_E / stops_F / stops_G
    // ================================
    const auto &stops = stops_D;
    return colormap_from_stops(stops, t);
}

static inline std::pair<float, float> auto_range_percentile(
    const std::vector<float> &v, float p_lo = 0.05f, float p_hi = 0.95f)
{
    std::vector<float> tmp;
    tmp.reserve(v.size());
    for(float x : v) {
        if(std::isfinite(x)) tmp.push_back(x);
    }
    if(tmp.empty()) return {0.0f, 1.0f};

    size_t lo_idx = static_cast<size_t>((tmp.size() - 1) * p_lo);
    size_t hi_idx = static_cast<size_t>((tmp.size() - 1) * p_hi);

    std::nth_element(tmp.begin(), tmp.begin() + lo_idx, tmp.end());
    float lo = tmp[lo_idx];
    std::nth_element(tmp.begin() + lo_idx + 1, tmp.begin() + hi_idx, tmp.end());
    float hi = tmp[hi_idx];

    if(!(std::isfinite(lo) && std::isfinite(hi)) || hi <= lo) {
        std::sort(tmp.begin(), tmp.end());
        lo = tmp.front();
        hi = tmp.back();
        if(hi <= lo) hi = lo + 1.0f;
    }
    return {lo, hi};
}

static inline bool write_png_rgb(const std::string &path, int w, int h,
                                const std::vector<uint8_t> &rgb)
{
    if(w <= 0 || h <= 0 || rgb.size() != static_cast<size_t>(w) * h * 3) return false;
    return stbi_write_png(path.c_str(), w, h, 3, rgb.data(), w * 3) != 0;
}

static inline bool write_tl_png(const std::string &path, int w, int h,
                                const std::vector<float> &tl_db)
{
    if(w <= 0 || h <= 0 || tl_db.size() != static_cast<size_t>(w) * h) return false;

    const auto [vmin, vmax] = auto_range_percentile(tl_db, 0.05f, 0.95f);

    std::vector<uint8_t> rgb(static_cast<size_t>(w) * h * 3);
    for(int iz = 0; iz < h; ++iz) {
        for(int ir = 0; ir < w; ++ir) {
            size_t idx = static_cast<size_t>(iz) * w + ir;
            float v = tl_db[idx];
            float t = (v - vmin) / (vmax - vmin);
            RGB c = colormap(t);
            size_t o = idx * 3;
            rgb[o+0] = c.r;
            rgb[o+1] = c.g;
            rgb[o+2] = c.b;
        }
    }

    return write_png_rgb(path, w, h, rgb);
}

} // namespace tl_png
