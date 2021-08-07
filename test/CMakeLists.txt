project(tests)

add_subdirectory(
        ${CMAKE_SOURCE_DIR}/external/googletest
        ${CMAKE_BINARY_DIR}/external/googletest
)

file(GLOB CPP_SOURCE_FILES ../src/*/*.cpp ../src/*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ../src/*.h ../src/*/*.h ../src/*/*/*.h)

file(GLOB CPP_TEST_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
# adding the Google_Tests_run target
add_executable(unittest ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES} ${CPP_TEST_SOURCE_FILES} core/util/utils.cpp)

set(GTEST_INCLUDE_DIR, "${CMAKE_BINARY_DIR}/external/googletest/googletest/include")
target_include_directories(unittest PUBLIC ${GTEST_INCLUDE_DIR})
target_link_libraries(unittest gtest gtest_main)

if (LAPACK_FOUND)
    target_link_libraries(unittest ${LAPACK_LIBRARIES})
endif ()


target_link_libraries(unittest ${MPI_LIBRARIES})
target_link_libraries(unittest yaml-cpp)

option(ENABLE_PROFILING "ENABLE_PROFILING" 0)
if (ENABLE_PROFILING)
    target_compile_options(unittest PRIVATE -pg)
    target_link_options(unittest PRIVATE -pg)
endif()

target_link_options(unittest PRIVATE -g -rdynamic)


if (ENABLE_HDF5)
    target_link_libraries(unittest ${_HDF5_LIBS})
    add_dependencies(unittest hdf5_local)
    target_link_libraries(unittest ${_HDF5_LIBS})
    target_include_directories(unittest PUBLIC
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/src
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/hl/src
            ${PROJECT_BINARY_DIR}/hdf5_local/build
    )
endif ()

#target_compile_definitions(unittest PRIVATE DEBUG=1)
#target_compile_options(unittest PRIVATE -O0)

include_directories("${CMAKE_SOURCE_DIR}")
include_directories($ENV{INCLUDE_PATH})
