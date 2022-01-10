cmake_minimum_required(VERSION 3.10)
project(M7)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif ()

add_subdirectory(${CMAKE_SOURCE_DIR}/external/spdlog-1.x)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
else ()
    message(FATAL_ERROR "Specified CMAKE_BUILD_TYPE is invalid")
endif()

find_package(LAPACK REQUIRED)


find_package(MPI REQUIRED COMPONENTS C)
message(STATUS "MPI_FOUND ${MPI_FOUND}")
if (MPI_FOUND)
    add_compile_definitions(ENABLE_MPI)
else ()
endif ()

#
# Download and build HDF5
#
set(HDF5_CMAKE_ARGS
        -DHDF5_ENABLE_PARALLEL=ON
        -DHDF5_BUILD_CPP_LIB=OFF
        -DHDF5_BUILD_JAVA=OFF
        -DHDF5_ENABLE_THREADSAFE=OFF
        -DCMAKE_C_FLAGS="-w") # we don't need to see compiler warnings from external dependencies

set(HDF5_USE_STATIC_LIBRARIES 1)
include(ExternalProject)
set(HDF5_VERSION "1.10.5")
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/hdf5_local/src)
ExternalProject_Add(hdf5_local
        URL https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-${HDF5_VERSION}.tar.gz
        URL_MD5 e115eeb66e944fa7814482415dd21cc4
        CMAKE_ARGS ${HDF5_CMAKE_ARGS}
        SOURCE_DIR ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${HDF5_VERSION}
        BINARY_DIR ${PROJECT_BINARY_DIR}/hdf5_local/build
        STAMP_DIR ${PROJECT_BINARY_DIR}/hdf5_local/stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/hdf5_local/tmp
        INSTALL_COMMAND ""
        )
include_directories(
        ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/src
        ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/hl/src
        ${PROJECT_BINARY_DIR}/hdf5_local/build
)
set(HDF5_LIBS
        ${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5.a
        ${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5_hl.a
        -ldl)

#
# Download and build static lib target of arpack-ng dependency of arpackpp
#
include(ExternalProject)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/arpack-ng/src)
ExternalProject_Add(arpack-ng
        URL https://github.com/opencollab/arpack-ng/archive/refs/heads/master.zip
        URL_MD5 88487b64fab4069f67e2b21c606d4602
        CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
        SOURCE_DIR ${PROJECT_BINARY_DIR}/arpack-ng/src/arpackpp
        BINARY_DIR ${PROJECT_BINARY_DIR}/arpack-ng/build
        STAMP_DIR ${PROJECT_BINARY_DIR}/arpack-ng/stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/arpack-ng/tmp
        INSTALL_COMMAND ""
        )

#
# Download and build arpack
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/arpackpp/src)
ExternalProject_Add(arpackpp
        URL https://github.com/robert-anderson/arpackpp/archive/refs/heads/master.zip
        URL_MD5 eb7eaaf9674e54b6c6b7785102a7c24b
        CMAKE_ARGS -DCMAKE_LIBRARY_PATH=${CMAKE_BINARY_DIR}/arpack-ng/build
        SOURCE_DIR ${PROJECT_BINARY_DIR}/arpackpp/src/arpackpp
        BINARY_DIR ${PROJECT_BINARY_DIR}/arpackpp/build
        STAMP_DIR ${PROJECT_BINARY_DIR}/arpackpp/stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/arpackpp/tmp
        INSTALL_COMMAND ""
        )

set(ARPACK_LIBS ${PROJECT_BINARY_DIR}/arpack-ng/build/libarpack.a)

#
# Command line args for M7
#
set(MBF_TYPE "0" CACHE STRING "Many-body basis function type index (0: determinants (default),
               1: determinant-permanent products, 2: permanents)")

if(${MBF_TYPE} LESS_EQUAL 2)
elseif (${MBF_TYPE} STREQUAL "fermion")
    set(MBF_TYPE 0)
elseif (${MBF_TYPE} STREQUAL "fermion-boson")
    set(MBF_TYPE 1)
elseif (${MBF_TYPE} STREQUAL "boson")
    set(MBF_TYPE 2)
else()
    message(FATAL_ERROR "Specified MBF_TYPE \"${MBF_TYPE}\" is invalid")
endif()

add_definitions(-DMBF_TYPE=${MBF_TYPE})

option(ENABLE_COMPLEX "ENABLE_COMPLEX" 0)
if (ENABLE_COMPLEX)
    add_compile_definitions(ENABLE_COMPLEX)
endif ()

option(ENABLE_ROUNDING_SPAWN_ATTEMPTS "ENABLE_ROUNDING_SPAWN_ATTEMPTS" 0)
if (ENABLE_ROUNDING_SPAWN_ATTEMPTS)
    add_compile_definitions(ENABLE_ROUNDING_SPAWN_ATTEMPTS)
endif ()

option(ENABLE_LOCAL_LOGGING "ENABLE_LOCAL_LOGGING" 0)
if (ENABLE_LOCAL_LOGGING)
    add_compile_definitions(ENABLE_LOCAL_LOGGING)
endif ()


include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")


include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/external/spdlog-1.x/include")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/yaml-cpp" yaml-cpp)
target_compile_options(yaml-cpp PRIVATE -w)
include_directories("${PROJECT_BINARY_DIR}/arpackpp/src/")

add_definitions(-DPROJECT_ROOT="${CMAKE_SOURCE_DIR}")
add_definitions(-Wall -Wextra -Werror)
add_definitions(-Wno-unused-function -Wno-unused-parameter)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(assets)
