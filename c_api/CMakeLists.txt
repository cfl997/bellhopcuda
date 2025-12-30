cmake_minimum_required(VERSION 3.15)

project(bhc_c_api LANGUAGES C CXX)

# 创建动态库
add_library(bhc_c_api SHARED
    src/bhc_c_api.cpp
    include/bhc/bhc_c_api.h
    include/bhc/bhc_tl.hpp
)

# 预处理器定义
target_compile_definitions(bhc_c_api PRIVATE 
    BHC_C_API_EXPORTS
    _CRT_SECURE_NO_WARNINGS # 消除 strncpy 等函数的安全警告
)

# 头文件搜索路径
target_include_directories(bhc_c_api
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/glm
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# 链接 bellhopcxx 静态库
# 假设 bellhopcxx 是静态库，如果不是，需要调整
if(TARGET bellhopcxx)
    # 获取 bellhopcxx 的依赖项
    get_target_property(bellhopcxx_type bellhopcxxstatic TYPE)
    if(bellhopcxx_type STREQUAL "STATIC_LIBRARY")
        target_link_libraries(bhc_c_api PRIVATE bellhopcxxstatic)
    else()
        message(FATAL_ERROR "bellhopcxxstatic 必须是静态库才能被链接到 bhc_c_api")
    endif()
elseif(TARGET bellhopcuda)
    # 如果 bellhopcuda 是静态库，直接链接
    get_target_property(bellhopcuda_type bellhopcuda TYPE)
    if(bellhopcuda_type STREQUAL "STATIC_LIBRARY")
        target_link_libraries(bhc_c_api PRIVATE bellhopcuda)
    else()
        message(FATAL_ERROR "bellhopcuda 必须是静态库才能被链接到 bhc_c_api")
    endif()
else()
    message(FATAL_ERROR "未找到 bellhopcxx 或 bellhopcuda 目标")
endif()

# C++ 标准
set_target_properties(bhc_c_api PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)

# Windows 下生成 .dll
if(WIN32)
    set_target_properties(bhc_c_api PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
    
    # 添加 UTF-8 支持
    target_compile_options(bhc_c_api PRIVATE /utf-8)
endif()

# 安装规则（可选）
install(TARGETS bhc_c_api
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/bhc
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)