include(${CMAKE_SOURCE_DIR}/cmake/build_or_find_hdf5.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/build_or_find_arpack.cmake)

add_library(M7_lib defs.h)
target_include_directories(M7_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

add_subdirectory(arnoldi)
add_subdirectory(basis)
add_subdirectory(bilinear)
add_subdirectory(caches)
add_subdirectory(config)
add_subdirectory(connection)
add_subdirectory(dynamics)
add_subdirectory(enumerator)
add_subdirectory(excitgen)
add_subdirectory(excitgen2)
add_subdirectory(excititer)
add_subdirectory(field)
add_subdirectory(foreach)
add_subdirectory(hamiltonian)
add_subdirectory(hash)
add_subdirectory(integrals)
add_subdirectory(io)
add_subdirectory(linalg)
add_subdirectory(mae)
add_subdirectory(nd)
add_subdirectory(observables)
add_subdirectory(parallel)
add_subdirectory(rdm)
add_subdirectory(sample)
add_subdirectory(sort)
add_subdirectory(subspace)
add_subdirectory(table)
add_subdirectory(util)
add_subdirectory(wavefunction)

add_subdirectory(
    "${CMAKE_SOURCE_DIR}/external/yaml-cpp"
    "${CMAKE_BINARY_DIR}/external/yaml-cpp"
)
target_compile_options(yaml-cpp PRIVATE -w)
target_link_libraries(M7_lib yaml-cpp)

add_subdirectory(
    ${CMAKE_SOURCE_DIR}/external/spdlog-1.x
    ${CMAKE_BINARY_DIR}/external/spdlog-1.x
)
target_link_libraries(M7_lib spdlog)

build_arpack(my_arpack)
target_link_libraries(M7_lib my_arpack)

find_package(LAPACK REQUIRED)
target_link_libraries(M7_lib LAPACK::LAPACK)

build_or_find_hdf5(M7_hdf5)
target_link_libraries(M7_lib M7_hdf5)

find_package(MPI REQUIRED COMPONENTS C)
target_compile_definitions(M7_lib PUBLIC ENABLE_MPI)
target_link_libraries(M7_lib MPI::MPI_C)

target_compile_definitions(M7_lib PUBLIC PROJECT_ROOT="${CMAKE_SOURCE_DIR}")


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

target_compile_definitions(M7_lib PUBLIC MBF_TYPE=${MBF_TYPE})

option(ENABLE_COMPLEX "ENABLE_COMPLEX")
if (ENABLE_COMPLEX)
    target_compile_definitions(M7_lib PRIVATE ENABLE_COMPLEX)
endif ()

option(ENABLE_LOCAL_LOGGING "ENABLE_LOCAL_LOGGING")
if (ENABLE_LOCAL_LOGGING)
    target_compile_definitions(M7_lib PRIVATE ENABLE_LOCAL_LOGGING)
endif ()

target_compile_options(M7_lib PUBLIC -Wno-cast-function-type)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(M7_lib -lgfortran)
    target_compile_options(M7_lib
        PUBLIC
             -mbmi -g
        PRIVATE
            -Wno-unused-function -Wno-unused-parameter
            -Wall -Wextra
    )
    if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
        target_compile_options(M7_lib
            PUBLIC
                -O3
        )
        target_compile_definitions(M7_lib
            PRIVATE
                NDEBUG
        )
    elseif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        target_compile_options(M7_lib
            PRIVATE
                -Werror
        )
    else ()
        message(FATAL_ERROR "Specified CMAKE_BUILD_TYPE is invalid")
    endif()
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} currently not supported.
                         Please contact the developers.")
endif()

get_target_property(compile_options M7_lib COMPILE_OPTIONS)
message(STATUS "The compile options are:")
message(STATUS "${compile_options}")
get_target_property(compile_definitions M7_lib COMPILE_DEFINITIONS)
message(STATUS "The compile definitions are:")
message(STATUS "${compile_definitions}")



set_source_files_properties(test.cpp
    PROPERTIES COMPILE_FLAGS
        -Wno-unused-function
)
