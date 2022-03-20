project(tests)

add_subdirectory(
        ${CMAKE_SOURCE_DIR}/external/googletest
        ${CMAKE_BINARY_DIR}/external/googletest
)

# adding the Google_Tests_run target
add_executable(unittest test.cpp)
target_include_directories(unittest
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
add_subdirectory(test_core)
target_link_libraries(unittest M7_lib)

target_link_libraries(unittest gtest gtest_main)

option(ENABLE_PROFILING "ENABLE_PROFILING" 0)
if (ENABLE_PROFILING)
    target_compile_options(unittest PRIVATE -pg)
    target_link_options(unittest PRIVATE -pg)
endif()

target_link_options(unittest PRIVATE -g -rdynamic)
