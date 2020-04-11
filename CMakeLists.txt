cmake_minimum_required(VERSION 3.10)
project(M7)
set(CMAKE_CXX_STANDARD 14)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

#SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
#IF( NOT EIGEN3_INCLUDE_DIR )
#    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
#ENDIF()
#INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

#INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/external/eigen")

add_subdirectory(${CMAKE_SOURCE_DIR}/external/CLI11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -fopenmp -march=native")

FIND_PACKAGE(MPI)
MESSAGE (STATUS "MPI_FOUND ${MPI_FOUND}")
if (MPI_FOUND)
    add_definitions(-DHAVE_MPI)
else()
endif()

include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/external/CLI11/include")

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(assets)
