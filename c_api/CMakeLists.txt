cmake_minimum_required(VERSION 3.15)

project(bhc_cpp_api LANGUAGES CXX)

add_library(bhc_cpp_api SHARED
    src/bhc_cpp_api.cpp
    include/bhc/bhc_cpp_api.hpp
)

target_compile_definitions(bhc_cpp_api PRIVATE
    _CRT_SECURE_NO_WARNINGS
    BHC_CPP_API_EXPORTS
)

target_include_directories(bhc_cpp_api
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/glm
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# 链接底层静态库
if(TARGET bellhopcxxstatic)
    target_link_libraries(bhc_cpp_api PRIVATE bellhopcxxstatic)
elseif(TARGET bellhopcudastatic)
    target_link_libraries(bhc_cpp_api PRIVATE bellhopcudastatic)
else()
    message(FATAL_ERROR "未找到 bellhopcxxstatic/bellhopcudastatic")
endif()

set_target_properties(bhc_cpp_api PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if(MSVC)
    target_compile_options(bhc_cpp_api PRIVATE /utf-8)
endif()
