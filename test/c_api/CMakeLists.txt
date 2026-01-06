cmake_minimum_required(VERSION 3.15)

project(bhc_cpp_api_tests LANGUAGES CXX)

add_executable(example_tl_cpp
    example_tl_cpp.cpp
)

add_executable(test_case_qvmt
    test_case_qvmt.cpp
)

add_executable(example_advanced
    example_advanced.cpp
)

# stb 是头文件库，无需额外链接

target_include_directories(example_tl_cpp PRIVATE
    ${CMAKE_SOURCE_DIR}/c_api/include
)

target_include_directories(test_case_qvmt PRIVATE
    ${CMAKE_SOURCE_DIR}/c_api/include
)

target_include_directories(example_advanced PRIVATE
    ${CMAKE_SOURCE_DIR}/c_api/include
)

target_link_libraries(example_tl_cpp PRIVATE
    bhc_cpp_api
)

target_link_libraries(test_case_qvmt PRIVATE
    bhc_cpp_api
)

target_link_libraries(example_advanced PRIVATE
    bhc_cpp_api
)

# ---- TEDSbellhop tests/examples ----
add_executable(example_teds_tl_cpp
    example_teds_tl_cpp.cpp
)

target_include_directories(example_teds_tl_cpp PRIVATE
    ${CMAKE_SOURCE_DIR}/TEDSbellhop/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(example_teds_tl_cpp PRIVATE
    TEDSbellhop
)

set_target_properties(example_teds_tl_cpp PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if(MSVC)
    target_compile_options(example_teds_tl_cpp PRIVATE /utf-8)

    add_custom_command(TARGET example_teds_tl_cpp POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:TEDSbellhop>
            $<TARGET_FILE_DIR:example_teds_tl_cpp>)
endif()

add_executable(test_teds_parallel
    test_teds_parallel.cpp
)

target_include_directories(test_teds_parallel PRIVATE
    ${CMAKE_SOURCE_DIR}/TEDSbellhop/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_teds_parallel PRIVATE
    TEDSbellhop
)

set_target_properties(test_teds_parallel PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if(MSVC)
    target_compile_options(test_teds_parallel PRIVATE /utf-8)

    add_custom_command(TARGET test_teds_parallel POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:TEDSbellhop>
            $<TARGET_FILE_DIR:test_teds_parallel>)
endif()

add_executable(test_teds_parallel_qssp
    test_teds_parallel_qssp.cpp
)

target_include_directories(test_teds_parallel_qssp PRIVATE
    ${CMAKE_SOURCE_DIR}/TEDSbellhop/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_teds_parallel_qssp PRIVATE
    TEDSbellhop
)

set_target_properties(test_teds_parallel_qssp PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if(MSVC)
    target_compile_options(test_teds_parallel_qssp PRIVATE /utf-8)

    add_custom_command(TARGET test_teds_parallel_qssp POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:TEDSbellhop>
            $<TARGET_FILE_DIR:test_teds_parallel_qssp>)
endif()


set_target_properties(example_tl_cpp PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

set_target_properties(test_case_qvmt PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

set_target_properties(example_advanced PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if(MSVC)
    target_compile_options(example_tl_cpp PRIVATE /utf-8)
    target_compile_options(test_case_qvmt PRIVATE /utf-8)
    target_compile_options(example_advanced PRIVATE /utf-8)

    add_custom_command(TARGET example_tl_cpp POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:bhc_cpp_api>
            $<TARGET_FILE_DIR:example_tl_cpp>)

    add_custom_command(TARGET test_case_qvmt POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:bhc_cpp_api>
            $<TARGET_FILE_DIR:test_case_qvmt>)

    add_custom_command(TARGET example_advanced POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:bhc_cpp_api>
            $<TARGET_FILE_DIR:example_advanced>)
endif()

