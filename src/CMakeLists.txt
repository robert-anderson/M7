project(src)

file(GLOB CPP_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ./*.h ./*/*.h ./*/*/*.h)

add_executable(M7 main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})

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

include_directories("${CMAKE_SOURCE_DIR}")
include_directories($ENV{INCLUDE_PATH})
