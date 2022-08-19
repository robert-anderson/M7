//
// Created by anderson on 27/06/2022.
//

#include "File.h"


bool hdf5::FileBase::is_hdf5(const str_t &fname) {
    return H5Fis_hdf5(fname.c_str()) > 0;
}

void hdf5::FileBase::require_is_hdf5(const str_t &fname) {
    REQUIRE_TRUE(is_hdf5(fname), "Specified file is not HDF5 format");
}