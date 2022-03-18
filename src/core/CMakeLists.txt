add_library(M7_lib defs.h)


target_include_directories(M7_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
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

target_link_libraries(M7_lib yaml-cpp)
target_link_libraries(M7_lib ${MPI_LIBRARIES})

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    target_link_libraries(M7_lib ${LAPACK_LIBRARIES})
endif()


target_link_libraries(M7_lib m7_hdf5)

target_link_libraries(M7_lib ${ARPACK_LIBS} -lgfortran)
add_dependencies(M7_lib arpack-ng)
add_dependencies(M7_lib arpackpp)
