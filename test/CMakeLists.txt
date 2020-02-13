project(tests)

add_subdirectory(
        ${CMAKE_SOURCE_DIR}/external/googletest
        ${CMAKE_BINARY_DIR}/external/googletest
)

file(GLOB CPP_SOURCE_FILES test.cpp ../src/*/*.cpp)
file(GLOB CPP_HEADER_FILES ../src/*/*.h)
file(GLOB CPP_TEST_SOURCE_FILES ./*/*.cpp)
# adding the Google_Tests_run target
add_executable(unittest ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES} ${CPP_TEST_SOURCE_FILES})

set(GTEST_INCLUDE_DIR, "${CMAKE_BINARY_DIR}/external/googletest/googletest/include")
target_include_directories(unittest PUBLIC ${GTEST_INCLUDE_DIR})
target_link_libraries(unittest gtest gtest_main)

target_compile_definitions(unittest PRIVATE DEBUG=1)
include_directories($ENV{INCLUDE_PATH})
SET(CMAKE_CXX_FLAGS "-fopenmp -march=native")