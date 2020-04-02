project(src)

file(GLOB CPP_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ./*.h ./*/*.h ./*/*/*.h)

add_executable(release main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})
add_executable(debug main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES} core/hash/HashMap.cpp core/hash/HashMap.h)

target_link_libraries (release CLI11::CLI11)
target_link_libraries (release ${MPI_LIBRARIES})
target_link_libraries (debug CLI11::CLI11)
target_link_libraries (debug ${MPI_LIBRARIES})

target_include_directories(release PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)
target_include_directories(debug PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)

target_compile_definitions(release PRIVATE NDEBUG=1)
target_compile_options(release PRIVATE -O2)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    target_link_libraries (release ${LAPACK_LIBRARIES})
    target_link_libraries (debug ${LAPACK_LIBRARIES})
endif()

include_directories($ENV{INCLUDE_PATH})
set(CMAKE_CXX_FLAGS "-fopenmp -march=native")