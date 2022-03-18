add_subdirectory(
        ${CMAKE_SOURCE_DIR}/external/googletest
        ${CMAKE_BINARY_DIR}/external/googletest
)

file(GLOB CPP_SOURCE_FILES ../src/*/*.cpp ../src/*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ../src/*.h ../src/*/*.h ../src/*/*/*.h)

file(GLOB CPP_TEST_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
file(GLOB CPP_TEST_HEADER_FILES ./*.h ./*/*.h ./*/*/*.h)
# adding the Google_Tests_run target
add_executable(unittest ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES}
        ${CPP_TEST_SOURCE_FILES} ${CPP_TEST_HEADER_FILES})

set(GTEST_INCLUDE_DIR, "${CMAKE_BINARY_DIR}/external/googletest/googletest/include")
target_include_directories(unittest PUBLIC ${GTEST_INCLUDE_DIR})
target_link_libraries(unittest gtest gtest_main)

if (LAPACK_FOUND)
    target_link_libraries(unittest ${LAPACK_LIBRARIES})
else()
    message(FATAL_ERROR "LAPACK implementation not found")
endif ()


target_link_libraries(unittest ${MPI_LIBRARIES})
target_link_libraries(unittest yaml-cpp)

option(ENABLE_PROFILING "ENABLE_PROFILING" 0)
if (ENABLE_PROFILING)
    target_compile_options(unittest PRIVATE -pg)
    target_link_options(unittest PRIVATE -pg)
endif()

target_link_options(unittest PRIVATE -g -rdynamic)

target_link_libraries(unittest m7_hdf5)

target_link_libraries(unittest ${ARPACK_LIBS} -lgfortran)
add_dependencies(unittest arpack-ng)
add_dependencies(unittest arpackpp)


include_directories("${CMAKE_SOURCE_DIR}")
include_directories($ENV{INCLUDE_PATH})
