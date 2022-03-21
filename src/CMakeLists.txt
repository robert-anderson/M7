add_executable(M7 main.cpp)
add_subdirectory(M7_lib)
target_link_libraries(M7 M7_lib)
set_target_properties(M7 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/bin"
)

