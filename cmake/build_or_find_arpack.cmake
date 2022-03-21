function(build_arpack target)
#
# Download and build static lib target of arpack-ng dependency of arpackpp
#
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/arpack-ng/src)
    ExternalProject_Add(arpack-ng
            URL https://github.com/robert-anderson/arpack-ng/archive/refs/heads/master.zip
            URL_MD5 88487b64fab4069f67e2b21c606d4602
            CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
            SOURCE_DIR ${PROJECT_BINARY_DIR}/arpack-ng/src/arpackpp
            BINARY_DIR ${PROJECT_BINARY_DIR}/arpack-ng/build
            STAMP_DIR ${PROJECT_BINARY_DIR}/arpack-ng/stamp
            TMP_DIR ${PROJECT_BINARY_DIR}/arpack-ng/tmp
            INSTALL_COMMAND ""
    )

#
# Download and build arpack
#
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/arpackpp/src)
    ExternalProject_Add(arpackpp
            URL https://github.com/robert-anderson/arpackpp/archive/refs/heads/master.zip
            URL_MD5 fc563ecfaf1ce07a1320a064ae3c02a8
            CMAKE_ARGS -DCMAKE_LIBRARY_PATH=${CMAKE_BINARY_DIR}/arpack-ng/build
            SOURCE_DIR ${PROJECT_BINARY_DIR}/arpackpp/src/arpackpp
            BINARY_DIR ${PROJECT_BINARY_DIR}/arpackpp/build
            STAMP_DIR ${PROJECT_BINARY_DIR}/arpackpp/stamp
            TMP_DIR ${PROJECT_BINARY_DIR}/arpackpp/tmp
            INSTALL_COMMAND ""
    )
    add_dependencies(arpackpp arpack-ng)

    add_library(${target} INTERFACE)
    add_dependencies(${target} arpackpp arpack-ng)
    target_include_directories(${target}
        INTERFACE
            ${PROJECT_BINARY_DIR}/arpackpp/src/arpackpp/include/
    )
    target_link_libraries(${target}
        INTERFACE
            ${PROJECT_BINARY_DIR}/arpack-ng/build/libarpack.a)
    target_link_options(${target}
        INTERFACE
            "-lgfortran"
    )
endfunction()

