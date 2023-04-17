cmake_minimum_required(VERSION 3.17)
project(M7 CXX C)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif ()

include(ExternalProject)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(assets)