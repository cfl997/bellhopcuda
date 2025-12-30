#ifndef BHC_C_API_H
#define BHC_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * 跨平台动态库导出/导入宏定义
 */
#if defined(_WIN32)
    #ifdef BHC_C_API_EXPORTS
        #define BHC_C_API __declspec(dllexport)
    #else
        #define BHC_C_API __declspec(dllimport)
    #endif
#else
    #define BHC_C_API __attribute__((visibility("default")))
#endif

/**
 * @brief 声速剖面(SSP)类型
 */
typedef enum {
    BHC_SSP_CLINEAR = 'C',  /**< C-线性插值 */
    BHC_SSP_N2LINEAR = 'N', /**< N²-线性插值 */
    BHC_SSP_CUBIC = 'S',    /**< 三次样条插值 */
    BHC_SSP_PCHIP = 'P'     /**< PCHIP插值 */
} BhcSSPType;

/**
 * @brief 海底地形插值类型
 */
typedef enum {
    BHC_BTY_LINEAR = 'L',   /**< 线性插值 */
    BHC_BTY_CURVATURE = 'C',/**< 曲率连续 */
    BHC_BTY_SPLINE = 'S'    /**< 样条插值 */
} BhcBathymetryType;

/**
 * @brief 波束类型
 */
typedef enum {
    BHC_BEAM_GEOM_HAT = 'G',     /**< 几何hat波束 */
    BHC_BEAM_GEOM_GAUSSIAN = 'g',/**< 几何高斯波束 */
    BHC_BEAM_CERVENY_CART = 'B', /**< Cerveny Cartesian波束 */
    BHC_BEAM_CERVENY_RAY = 'b'   /**< Cerveny 射线中心波束 */
} BhcBeamType;

/**
 * @brief 边界类型
 */
typedef enum {
    BHC_BOUNDARY_VACUUM = 'V',   /**< 真空边界 */
    BHC_BOUNDARY_RIGID = 'R',    /**< 刚性边界 */
    BHC_BOUNDARY_ACOUSTIC = 'A', /**< 声学半空间 */
    BHC_BOUNDARY_FILE = 'F'      /**< 文件定义边界 */
} BhcBoundaryType;

/**
 * @brief 运行类型
 */
typedef enum {
    BHC_RUN_COHERENT = 'C',     /**< 相干TL计算 */
    BHC_RUN_INCOHERENT = 'I',   /**< 非相干TL计算 */
    BHC_RUN_SEMI_COHERENT = 'S',/**< 半相干TL计算 */
    BHC_RUN_EIGENRAYS = 'E',    /**< 特征声线计算 */
    BHC_RUN_ARRIVALS = 'A',     /**< 本征声线(ASCII输出) */
    BHC_RUN_ARRIVALS_BIN = 'a'  /**< 本征声线(二进制输出) */
} BhcRunType;

/**
 * @brief 计算维度
 */
typedef enum {
    BHC_DIM_2D = 2,   /**< 2D计算 */
    BHC_DIM_3D = 3,   /**< 3D计算 */
    BHC_DIM_NX2D = 4  /**< Nx2D计算 */
} BhcDimension;

/**
 * @struct BhcSSP
 * @brief 定义声速剖面 (Sound Speed Profile)
 */
typedef struct {
    int n_points;      /**< 剖面点数 */
    double* z;         /**< 深度数组 (m), 大小 `n_points`，必须按深度递增排序 */
    double* c;         /**< 声速数组 (m/s), 大小 `n_points` */
    double* rho;       /**< 密度数组 (g/cm³), 大小 `n_points` */
    double* alpha;     /**< 声吸收系数 (dB/λ), 大小 `n_points`，可选，可为NULL */
    BhcSSPType type;   /**< SSP类型 */
} BhcSSP;

/**
 * @struct BhcBathymetry
 * @brief 定义海底地形
 */
typedef struct {
    int n_points;           /**< 地形点数 */
    double* r;              /**< 水平距离数组 (m), 大小 `n_points` */
    double* z;              /**< 深度数组 (m), 大小 `n_points` */
    BhcBathymetryType type; /**< 地形插值类型 */
} BhcBathymetry;

/**
 * @struct BhcBeamOptions
 * @brief 波束追踪选项
 */
typedef struct {
    int n_beams;          /**< 波束数量 */
    double angle_min;     /**< 最小发射角 (度) */
    double angle_max;     /**< 最大发射角 (度) */
    BhcBeamType beam_type;/**< 波束类型 */
    double step_size;     /**< 步长 (m), 如果<=0则自动计算 */
    double max_range;     /**< 最大计算距离 (m) */
    double beam_width;    /**< 波束宽度 (度), 用于波束窗口 */
    int max_reflections;  /**< 最大反射次数 */
    int max_refractions;  /**< 最大折射次数 */
} BhcBeamOptions;

/**
 * @struct BhcBoundary
 * @brief 边界条件设置
 */
typedef struct {
    BhcBoundaryType surface_type;  /**< 海面类型 */
    double surface_roughness;      /**< 海面粗糙度 (m), 0表示平滑 */
    BhcBoundaryType bottom_type;   /**< 海底类型 */
    double bottom_rho;             /**< 海底密度 (g/cm³) */
    double bottom_cp;              /**< 海底压缩波速 (m/s) */
    double bottom_cs;              /**< 海底剪切波速 (m/s), 0表示无剪切波 */
    double bottom_ap;              /**< 海底压缩波衰减 (dB/λ) */
    double bottom_as;              /**< 海底剪切波衰减 (dB/λ) */
    double bottom_roughness;       /**< 海底粗糙度 (m), 0表示平滑 */
} BhcBoundary;

/**
 * @struct BhcConfig
 * @brief TL计算的整体配置
 */
typedef struct {
    // 基本设置
    const char* title;     /**< 计算标题，将写入输出文件 */
    BhcDimension dimension;/**< 计算维度 */
    
    // 运行类型
    BhcRunType run_type;   /**< 运行类型 */
    
    // 源设置
    double src_depth;      /**< 源深度 (m) */
    double src_range;      /**< 源水平距离 (m), 通常为0 */
    double src_bearing;    /**< 源方位角 (度), 用于3D */
    
    // 接收器设置
    int n_ranges;          /**< 接收器距离点数 */
    int n_depths;          /**< 接收器深度点数 */
    int n_bearings;        /**< 接收器方位角数 (3D) */
    double range_min;      /**< 最小接收距离 (m) */
    double range_max;      /**< 最大接收距离 (m) */
    double depth_min;      /**< 最小接收深度 (m) */
    double depth_max;      /**< 最大接收深度 (m) */
    double bearing_min;    /**< 最小方位角 (度), 用于3D */
    double bearing_max;    /**< 最大方位角 (度), 用于3D */
    
    // 频率设置
    double frequency;      /**< 信号频率 (Hz) */
    int n_freqs;           /**< 频率点数 (用于宽带计算) */
    double* freqs;         /**< 频率数组 (Hz), 大小 `n_freqs`，用于宽带计算 */
    
    // 波束选项
    BhcBeamOptions beam;   /**< 波束追踪选项 */
    
    // 边界条件
    BhcBoundary boundary;  /**< 边界条件设置 */
    
    // 输出选项
    int output_rays;       /**< 是否输出射线轨迹: 0=否, 1=是 */
    int output_eigenrays;  /**< 是否输出本征声线: 0=否, 1=是 */
    int output_arrivals;   /**< 是否输出到达结构: 0=否, 1=是 */
    
    // 高级选项
    double eps_multiplier; /**< 步长乘数，影响计算精度和速度 */
    int max_steps;         /**< 最大步数，0表示使用默认值 */
    
} BhcConfig;

/**
 * @struct BhcTLResult
 * @brief 传输损失(TL)计算结果
 */
typedef struct {
    int n_freqs;        /**< 频率点数 */
    int n_ranges;       /**< 距离点数 */
    int n_depths;       /**< 深度点数 */
    int n_bearings;     /**< 方位角数 (3D) */
    
    double* freqs;      /**< 频率数组 (Hz), 大小 `n_freqs` */
    double* ranges;     /**< 距离数组 (m), 大小 `n_ranges` */
    double* depths;     /**< 深度数组 (m), 大小 `n_depths` */
    double* bearings;   /**< 方位角数组 (度), 大小 `n_bearings` (3D) */
    
    /**
     * TL 值数组 (dB)
     * 内存布局: [freq][bearing][range][depth]
     * 访问方式: ((freq * n_bearings + bearing) * n_ranges + range) * n_depths + depth
     */
    float* tl;          
    
    // 内部使用，用户不应直接访问
    void* _internal;    
} BhcTLResult;

/**
 * @brief 计算传输损失(TL)
 * 
 * 这是主要的计算函数，它接受配置、SSP和地形数据，执行计算并返回结果。
 * 返回的指针必须使用 bhc_free_tl_result 释放。
 *
 * @param config 计算配置
 * @param ssp 声速剖面
 * @param bty 海底地形 (可为NULL，表示平坦海底)
 * @return 指向 BhcTLResult 结构体的指针，失败时返回 NULL
 */
BHC_C_API BhcTLResult* bhc_compute_tl(const BhcConfig* config, 
                                     const BhcSSP* ssp, 
                                     const BhcBathymetry* bty);

/**
 * @brief 释放 BhcTLResult 占用的内存
 * 
 * @param result 由 bhc_compute_tl 返回的结果指针
 */
BHC_C_API void bhc_free_tl_result(BhcTLResult* result);

/**
 * @brief 获取最后一条错误信息
 * 
 * @return 错误信息字符串，如果没有错误则返回 NULL
 */
BHC_C_API const char* bhc_get_last_error(void);

/**
 * @brief 获取库版本信息
 * 
 * @return 版本信息字符串
 */
BHC_C_API const char* bhc_get_version(void);

/**
 * @brief 初始化 Bellhop 环境
 * 
 * 在调用其他 API 函数前，应该先调用此函数进行初始化。
 * 如果成功返回 0，失败返回非零错误码。
 */
BHC_C_API int bhc_initialize(void);

/**
 * @brief 清理 Bellhop 环境
 * 
 * 在程序退出前调用此函数释放所有资源。
 */
BHC_C_API void bhc_cleanup(void);

#ifdef __cplusplus
} // extern "C"

// C++11 封装头文件
#if __cplusplus >= 201103L
#include "bhc/bhc_tl.hpp"
#endif // C++11

#endif // __cplusplus

#endif // BHC_C_API_H