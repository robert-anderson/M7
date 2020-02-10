project(src)

file(GLOB CPP_SOURCE_FILES ./*/*.cpp)
file(GLOB CPP_HEADER_FILES ./*/*.h)

add_executable(release main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})
add_executable(debug main.cpp ${CPP_HEADER_FILES} ${CPP_SOURCE_FILES})

target_include_directories(release PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)
target_include_directories(debug PRIVATE ${CMAKE_SOURCE_DIR}/external/CLI11/include)

target_compile_definitions(debug PRIVATE DEBUG=1)
include_directories($ENV{INCLUDE_PATH})
SET(CMAKE_CXX_FLAGS "-fopenmp -march=native")