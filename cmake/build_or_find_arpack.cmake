function(add_arpack target)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/arpack-ng/src)
    ExternalProject_Add(arpack-ng
        URL https://github.com/robert-anderson/arpack-ng/archive/refs/heads/master.zip
        URL_MD5 88487b64fab4069f67e2b21c606d4602
        CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
        SOURCE_DIR ${PROJECT_BINARY_DIR}/arpack-ng/src/arpack-ng
        BINARY_DIR ${PROJECT_BINARY_DIR}/arpack-ng/build
        STAMP_DIR ${PROJECT_BINARY_DIR}/arpack-ng/stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/arpack-ng/tmp
        INSTALL_COMMAND ""
    )
add_dependencies(${target} arpack-ng)
#target_link_options(arpack-ng INTERFACE "-lgfortran")
target_include_directories(${target} PUBLIC ${CMAKE_SOURCE_DIR}/external/arpackpp/include)
target_link_libraries(${target} ${PROJECT_BINARY_DIR}/arpack-ng/build/libarpack.a)
endfunction()