function(build_or_find_hdf5 target)
    option(BUILD_HDF5 OFF)
    if(BUILD_HDF5)
#
# Download and build HDF5
#
        set(OWN_HDF5_ROOT ${PROJECT_BINARY_DIR}/lib/hdf5)
        set(OWN_HDF5_BUILD ${PROJECT_BINARY_DIR}/hdf5_local/)
        set(HDF5_USE_STATIC_LIBRARIES 1)
        set(HDF5_VERSION "1.10.5")

        file(MAKE_DIRECTORY ${OWN_HDF5_ROOT})
        file(MAKE_DIRECTORY ${OWN_HDF5_BUILD})
        ExternalProject_Add(built_hdf5
                URL
                    https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-${HDF5_VERSION}.tar.gz
                URL_MD5
                    e115eeb66e944fa7814482415dd21cc4
                CMAKE_ARGS
                    -DHDF5_ENABLE_PARALLEL=ON
                    -DHDF5_BUILD_CPP_LIB=OFF
                    -DHDF5_BUILD_JAVA=OFF
                    -DCMAKE_INSTALL_PREFIX=${OWN_HDF5_ROOT}
                    -DHDF5_ENABLE_THREADSAFE=OFF
                    -DHDF5_BUILD_TOOLS=OFF
                    -DHDF5_BUILD_EXAMPLES=OFF
                    -DBUILD_TESTING=OFF
                    -DCMAKE_C_FLAGS="-w"
                SOURCE_DIR
                    ${OWN_HDF5_BUILD}/src/hdf5-${HDF5_VERSION}
                BINARY_DIR
                    ${OWN_HDF5_BUILD}/build
                STAMP_DIR
                    ${OWN_HDF5_BUILD}/stamp
                TMP_DIR
                    ${OWN_HDF5_BUILD}/tmp
                INSTALL_DIR
                    ""
                INSTALL_COMMAND
                    make install
        )
        add_library(${target} INTERFACE)
        add_dependencies(${target} built_hdf5)
        target_include_directories(${target}
            INTERFACE
                ${OWN_HDF5_ROOT}/include
        )
        target_link_libraries(${target}
            INTERFACE
                ${OWN_HDF5_ROOT}/lib/libhdf5_hl.a
                ${OWN_HDF5_ROOT}/lib/libhdf5.a
        )
    else()
        find_package(HDF5 COMPONENTS C HL)
        if (NOT ${HDF5_FOUND})
            message(FATAL_ERROR "HDF5 not found. Either specify the environment "
                "variable HDF5_ROOT to guide to a valid HDF5 installation, "
                "or use the option `-DBUILD_HDF5=ON` to download "
                "and compile HDF5 automatically."
            )
        endif()

        add_library(${target} INTERFACE)
        target_link_libraries(${target} INTERFACE hdf5::hdf5 hdf5::hdf5_hl)
    endif()
    target_link_libraries(${target} INTERFACE ${CMAKE_DL_LIBS})
endfunction()
