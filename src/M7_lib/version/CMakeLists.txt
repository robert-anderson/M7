# we could just set the version as a compile definition but this would require recompiling all files once the commit
# hash changes, which is inefficient and promotes infrequent committing. instead a separate target is created for
# this purpose.
# also, this special target is used to bring in the correct compilation timestamp. the __DATE__, __TIME__ macros defined
# by the compiler will only change when the files that call them are recompiled, so cmake variables are a better solution

add_library(M7_version STATIC version.h)

# adapted from https://www.mattkeeter.com/blog/2018-01-06-versioning/
execute_process(COMMAND git log --pretty=format:%h -n 1 OUTPUT_VARIABLE GIT_REV ERROR_QUIET)

string(TIMESTAMP COMPILATION_TIMESTAMP "%Y-%m-%d %H:%M:%S")

# Check whether we got any revision (which isn't
# always the case, e.g. when someone downloaded a zip
# file from Github instead of a checkout)
if ("${GIT_REV}" STREQUAL "")
    set(GIT_REV "N/A")
endif()

target_compile_definitions(M7_version PRIVATE VERSION="${GIT_REV}")
target_compile_definitions(M7_version PRIVATE COMPILATION_TIMESTAMP="${COMPILATION_TIMESTAMP}")
target_sources(M7_version PRIVATE version.cpp version.h)

set_target_properties(M7_version PROPERTIES LINKER_LANGUAGE C)
add_custom_command(TARGET M7_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt)
