project(src)

file(GLOB CPP_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ./*.h ./*/*.h ./*/*/*.h)

add_executable(release main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})
add_executable(debug main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})
add_executable(prof main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})

target_link_libraries (release yaml-cpp)
target_link_libraries (release ${MPI_LIBRARIES})
target_link_libraries (debug yaml-cpp)
target_link_libraries (debug ${MPI_LIBRARIES})
target_link_libraries (prof yaml-cpp)
target_link_libraries (prof ${MPI_LIBRARIES})

target_compile_definitions(release PRIVATE NDEBUG=1)
target_compile_options(release PRIVATE -O3)
target_link_options(release PRIVATE -g -rdynamic)

target_compile_definitions(debug PRIVATE VERBOSE=1)
target_link_options(debug PRIVATE)
target_link_options(debug PRIVATE -g -rdynamic)

target_compile_definitions(prof PRIVATE NDEBUG=1)
target_compile_options(prof PRIVATE -pg)
target_link_options(prof PRIVATE -pg)
target_link_options(prof PRIVATE -g -rdynamic)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    target_link_libraries (release ${LAPACK_LIBRARIES})
    target_link_libraries (debug ${LAPACK_LIBRARIES})
    target_link_libraries (prof ${LAPACK_LIBRARIES})
endif()

if (ENABLE_HDF5)
    target_link_libraries(release ${_HDF5_LIBS})
    target_link_libraries(debug ${_HDF5_LIBS})
    target_link_libraries(prof ${_HDF5_LIBS})
    add_dependencies(release hdf5_local)
    add_dependencies(debug hdf5_local)
    add_dependencies(prof hdf5_local)
    target_link_libraries(release ${_HDF5_LIBS})
    target_link_libraries(debug ${_HDF5_LIBS})
    target_link_libraries(prof ${_HDF5_LIBS})
    target_include_directories(release PUBLIC
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/src
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/hl/src
            ${PROJECT_BINARY_DIR}/hdf5_local/build
            )
    target_include_directories(debug PUBLIC
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/src
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/hl/src
            ${PROJECT_BINARY_DIR}/hdf5_local/build
            )
    target_include_directories(prof PUBLIC
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/src
            ${PROJECT_BINARY_DIR}/hdf5_local/src/hdf5-${_HDF5_VERSION}/hl/src
            ${PROJECT_BINARY_DIR}/hdf5_local/build
            )
endif ()

include_directories("${CMAKE_SOURCE_DIR}")
include_directories($ENV{INCLUDE_PATH})
