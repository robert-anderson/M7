project(tests)

add_subdirectory(
        ${CMAKE_SOURCE_DIR}/external/googletest
        ${CMAKE_BINARY_DIR}/external/googletest
)

file(GLOB CPP_SOURCE_FILES ../src/*/*.cpp ../src/*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ../src/*.h ../src/*/*.h ../src/*/*/*.h)

file(GLOB CPP_TEST_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
# adding the Google_Tests_run target
add_executable(unittest ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES} ${CPP_TEST_SOURCE_FILES})

set(GTEST_INCLUDE_DIR, "${CMAKE_BINARY_DIR}/external/googletest/googletest/include")
target_include_directories(unittest PUBLIC ${GTEST_INCLUDE_DIR})
target_link_libraries(unittest gtest gtest_main)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    target_link_libraries (unittest ${LAPACK_LIBRARIES})
endif()
target_link_libraries (unittest ${MPI_LIBRARIES})

FIND_PACKAGE(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
    message(${HDF5_INCLUDE_DIR})
    include_directories(${HDF5_INCLUDE_DIR})
    link_directories(${HDF5_LIBRARY_DIRS})
    target_link_libraries (unittest ${HDF5_CXX_LIBRARIES})
    if (${HDF5_IS_PARALLEL})
        target_compile_definitions(unittest PRIVATE HAVE_HDF5_PARALLEL)
    endif ()
endif ()

target_compile_definitions(unittest PRIVATE DEBUG=1)
#target_compile_definitions(unittest PRIVATE VERBOSE=1)
target_compile_options(unittest PRIVATE -O0)

include_directories($ENV{INCLUDE_PATH})
SET(CMAKE_CXX_FLAGS "-fopenmp -march=native")