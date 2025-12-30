#pragma once

#include <vector>
#include <memory>
#include <stdexcept>
#include "bhc_c_api.h"

namespace bhc {

/**
 * @brief 表示 TL 计算结果的 C++11 封装类
 * 
 * 提供类似多维数组的访问接口，支持 2D/3D 和频率分层的 TL 数据。
 * 内存布局：最内层是 depth，然后是 range，bearing，最后是频率。
 */
class TLResult {
public:
    // 禁用复制构造和赋值
    TLResult(const TLResult&) = delete;
    TLResult& operator=(const TLResult&) = delete;

    // 移动构造和赋值
    TLResult(TLResult&&) = default;
    TLResult& operator=(TLResult&&) = default;

    ~TLResult() = default;

    /**
     * @brief 获取 TL 值 (2D 情况)
     * @param range_idx 距离索引 [0, width()-1]
     * @param depth_idx 深度索引 [0, height()-1]
     * @return TL 值 (dB)
     */
    float at(int range_idx, int depth_idx) const {
        return at(0, 0, range_idx, depth_idx);
    }

    /**
     * @brief 获取 TL 值 (完整维度)
     * @param freq_idx 频率索引 [0, layers()-1]
     * @param bearing_idx 方位角索引 [0, faces()-1]
     * @param range_idx 距离索引 [0, width()-1]
     * @param depth_idx 深度索引 [0, height()-1]
     * @return TL 值 (dB)
     */
    float at(int freq_idx, int bearing_idx, int range_idx, int depth_idx) const {
        if (!data_ || !check_bounds(freq_idx, bearing_idx, range_idx, depth_idx)) {
            throw std::out_of_range("TLResult: 索引越界");
        }
        return data_[calc_index(freq_idx, bearing_idx, range_idx, depth_idx)];
    }

    // 维度信息
    int width() const noexcept { return width_; }     // 距离方向点数
    int height() const noexcept { return height_; }   // 深度方向点数
    int faces() const noexcept { return faces_; }     // 方位角数 (2D=1)
    int layers() const noexcept { return layers_; }   // 频率数

    // 获取原始数据指针 (谨慎使用)
    const float* data() const noexcept { return data_.get(); }
    float* data() noexcept { return data_.get(); }

    // 从 C 接口结果创建 (转移所有权)
    static TLResult from_c_result(BhcTLResult* c_result) {
        if (!c_result) {
            throw std::invalid_argument("TLResult: 输入指针不能为空");
        }
        return TLResult(c_result);
    }

private:
    // 私有构造函数，只能通过 from_c_result 创建
    explicit TLResult(BhcTLResult* c_result) 
        : width_(c_result->n_ranges),
          height_(c_result->n_depths),
          faces_(c_result->n_bearings > 0 ? c_result->n_bearings : 1),
          layers_(c_result->n_freqs > 0 ? c_result->n_freqs : 1),
          data_(c_result->tl, [](float* p) { 
              // 自定义删除器，不实际删除内存，由 C API 管理
          }) {
        // 转移所有权后清空 C 结构体
        c_result->tl = nullptr;
        bhc_free_tl_result(c_result);
    }

    // 计算一维索引
    size_t calc_index(int freq_idx, int bearing_idx, int range_idx, int depth_idx) const noexcept {
        return ((freq_idx * faces_ + bearing_idx) * width_ + range_idx) * height_ + depth_idx;
    }

    // 检查索引是否越界
    bool check_bounds(int freq_idx, int bearing_idx, int range_idx, int depth_idx) const noexcept {
        return (freq_idx >= 0 && freq_idx < layers_ &&
                bearing_idx >= 0 && bearing_idx < faces_ &&
                range_idx >= 0 && range_idx < width_ &&
                depth_idx >= 0 && depth_idx < height_);
    }

    int width_;     // 距离方向点数
    int height_;    // 深度方向点数
    int faces_;     // 方位角数 (2D=1)
    int layers_;    // 频率数
    std::shared_ptr<float> data_;  // 共享指针管理内存
};

/**
 * @brief 计算 TL 的 C++11 封装函数
 * 
 * @param config 配置参数
 * @param ssp 声速剖面
 * @param bty 海底地形 (可选)
 * @return TLResult TL 计算结果
 * @throw std::runtime_error 计算失败时抛出异常
 */
inline TLResult compute_tl(const BhcConfig& config, const BhcSSP& ssp, const BhcBathymetry* bty = nullptr) {
    // 调用 C API
    BhcTLResult* result = bhc_compute_tl(&config, &ssp, bty);
    
    // 检查错误
    if (!result) {
        const char* err = bhc_get_last_error();
        throw std::runtime_error(err ? err : "未知错误: bhc_compute_tl 返回空指针");
    }
    
    // 转换为 C++ 对象并返回
    return TLResult::from_c_result(result);
}

} // namespace bhc
