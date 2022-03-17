add_executable(M7 main.cpp main.cpp defs.h)
add_subdirectory(core)

target_link_libraries (M7 yaml-cpp)
target_link_libraries (M7 ${MPI_LIBRARIES})

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    target_link_libraries (M7 ${LAPACK_LIBRARIES})
endif()


target_link_libraries(M7 m7_hdf5)

target_link_libraries(M7 ${ARPACK_LIBS} -lgfortran)
add_dependencies(M7 arpack-ng)
add_dependencies(M7 arpackpp)
