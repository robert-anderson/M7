add_executable(M7 main.cpp)
add_subdirectory(M7_lib)
target_link_libraries(M7 M7_lib)

