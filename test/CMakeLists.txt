project(tests)

add_subdirectory(
        ${CMAKE_SOURCE_DIR}/external/googletest
        ${CMAKE_BINARY_DIR}/external/googletest
)

# adding the Google_Tests_run target
add_executable(unittest test.cpp)
target_include_directories(unittest
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
add_subdirectory(test_core)
target_link_libraries(unittest M7_lib)

target_link_libraries(unittest gtest gtest_main)

option(ENABLE_PROFILING "ENABLE_PROFILING")
if (ENABLE_PROFILING)
    target_compile_options(unittest PRIVATE -pg)
    target_link_options(unittest PRIVATE -pg)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(unittest
        PUBLIC
            -g
        PRIVATE
            -Wno-unused-function -Wno-unused-parameter
            -rdynamic
            -Wall -Wextra
    )
    if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
        target_compile_options(unittest
        )
    elseif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        target_compile_options(unittest
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
