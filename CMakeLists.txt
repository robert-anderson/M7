cmake_minimum_required(VERSION 3.10)
project(M7)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/external/spdlog-1.x)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

find_package(LAPACK REQUIRED)


find_package(MPI REQUIRED COMPONENTS C)
message(STATUS "MPI_FOUND ${MPI_FOUND}")
if (MPI_FOUND)
    add_compile_definitions(ENABLE_MPI)
else()
endif()

set(ExternalProjectCMakeArgs
        -DHDF5_ENABLE_PARALLEL=ON
        -DHDF5_BUILD_CPP_LIB=OFF
        -DHDF5_BUILD_JAVA=OFF
        -DHDF5_ENABLE_THREADSAFE=OFF
        -DCMAKE_C_FLAGS="-Wno-implicit-function-declaration")

option(FORCE_HDF5_INSTALL "Force the download and installation of HDF5" 1)
option(ENABLE_HDF5 "Enable HDF5 support" 1)
if(ENABLE_HDF5)
    set(HDF5_USE_STATIC_LIBRARIES 1)
    #find_package(HDF5)
    #if(HDF5_FOUND)
    #    include_directories(${HDF5_INCLUDE_DIR})
    #    set(_HDF5_LIBS hdf5)
    #else()
        # we did not find it so we fetch it from the web, configure, compile, and link
        include(ExternalProject)
        set(_HDF5_VERSION "1.10.5")
        file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/hdf5_local/src)
        ExternalProject_Add(hdf5_local
                URL https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-${_HDF5_VERSION}.tar.gz
                #URL /home/rja/Downloads/hdf5-${_HDF5_VERSION}.tar.gz
                URL_MD5 e115eeb66e944fa7814482415dd21cc4
                CMAKE_ARGS ${ExternalProjectCMakeArgs}
                SOURCE_DIR ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}
                BINARY_DIR ${PROJECT_BINARY_DIR}/hdf5_local/build
                STAMP_DIR  ${PROJECT_BINARY_DIR}/hdf5_local/stamp
                TMP_DIR    ${PROJECT_BINARY_DIR}/hdf5_local/tmp
                INSTALL_COMMAND ""
                )
        include_directories(
                ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/src
                ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/hl/src
                ${PROJECT_BINARY_DIR}/hdf5_local/build
        )
        set(_HDF5_LIBS
                ${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5.a
                ${PROJECT_BINARY_DIR}/hdf5_local/build/bin/libhdf5_hl.a -ldl)
    #endif()
endif()

#
# Command line args for M7
#
option(ENABLE_BOSONS "ENABLE_BOSONS" 0)
if (ENABLE_BOSONS)
    add_compile_definitions(ENABLE_BOSONS)
endif()

option(ENABLE_COMPLEX "ENABLE_COMPLEX" 0)
if (ENABLE_COMPLEX)
    add_compile_definitions(ENABLE_COMPLEX)
endif()

option(ENABLE_MEVS "ENABLE_MEVS" 0)
if (ENABLE_MEVS)
    add_compile_definitions(ENABLE_MEVS)
endif()

option(ENABLE_ROUNDING_SPAWN_ATTEMPTS "ENABLE_ROUNDING_SPAWN_ATTEMPTS" 0)
if (ENABLE_ROUNDING_SPAWN_ATTEMPTS)
    add_compile_definitions(ENABLE_ROUNDING_SPAWN_ATTEMPTS)
endif()

option(ENABLE_LOCAL_LOGGING "ENABLE_LOCAL_LOGGING" 0)
if (ENABLE_LOCAL_LOGGING)
    add_compile_definitions(ENABLE_LOCAL_LOGGING)
endif()


include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")


include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/external/spdlog-1.x/include")
include_directories("${CMAKE_SOURCE_DIR}/external/lambda-lanczos/include")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/yaml-cpp" yaml-cpp)

add_definitions(-DPROJECT_ROOT="${CMAKE_SOURCE_DIR}")
add_definitions(-Wall -Wextra -Werror)
add_definitions(-Wno-unused-function -Wno-unused-parameter -Wno-unused-value)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(assets)
