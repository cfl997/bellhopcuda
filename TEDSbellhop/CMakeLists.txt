cmake_minimum_required(VERSION 3.15)

project(TEDSbellhop LANGUAGES CXX)

add_library(TEDSbellhop SHARED
    src/TEDSbellhop.cpp
    include/TEDSbellhop/TEDSbellhop.hpp
)

target_compile_definitions(TEDSbellhop PRIVATE
    _CRT_SECURE_NO_WARNINGS
    TEDSBELLHOP_EXPORTS
)

target_include_directories(TEDSbellhop
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/glm
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# 链接底层静态库（与 c_api 保持一致）
if(TARGET bellhopcxxstatic)
    target_link_libraries(TEDSbellhop PRIVATE bellhopcxxstatic)
elseif(TARGET bellhopcudastatic)
    target_link_libraries(TEDSbellhop PRIVATE bellhopcudastatic)
else()
    message(FATAL_ERROR "未找到 bellhopcxxstatic/bellhopcudastatic")
endif()

set_target_properties(TEDSbellhop PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if(MSVC)
    target_compile_options(TEDSbellhop PRIVATE /utf-8)
endif()

