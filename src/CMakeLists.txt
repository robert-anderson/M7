project(src)

file(GLOB CPP_SOURCE_FILES ./*.cpp ./*/*.cpp ./*/*/*.cpp)
file(GLOB CPP_HEADER_FILES ./*.h ./*/*.h ./*/*/*.h)

add_executable(release main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})
add_executable(debug main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES} core/data/Buffer.cpp core/data/Buffer.h core/data/BufferWindow.cpp core/data/BufferWindow.h core/data/Table.cpp core/data/Table.h core/data/FieldBase.cpp core/data/FieldBase.h core/data/Field.cpp core/data/Field.h core/nd/NdArrayFormat.h core/data/BitsetField.cpp core/data/BitsetField.h core/data/NumericField.h)
add_executable(prof main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})

target_link_libraries (release CLI11::CLI11)
target_link_libraries (release ${MPI_LIBRARIES})
target_link_libraries (debug CLI11::CLI11)
target_link_libraries (debug ${MPI_LIBRARIES})
target_link_libraries (prof CLI11::CLI11)
target_link_libraries (prof ${MPI_LIBRARIES})

target_include_directories(release PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)
target_include_directories(debug PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)
target_include_directories(prof PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)

target_compile_definitions(release PRIVATE NDEBUG=1)
target_compile_options(release PRIVATE -O3)

target_compile_definitions(debug PRIVATE VERBOSE=1)

target_compile_definitions(prof PRIVATE NDEBUG=1)
target_compile_options(prof PRIVATE -O2)
target_compile_options(prof PRIVATE -pg)
target_link_options(prof PRIVATE -pg)

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    target_link_libraries (release ${LAPACK_LIBRARIES})
    target_link_libraries (debug ${LAPACK_LIBRARIES})
    target_link_libraries (prof ${LAPACK_LIBRARIES})
endif()

include_directories("${CMAKE_SOURCE_DIR}")
include_directories($ENV{INCLUDE_PATH})
