//
// Created by anderson on 27/06/2022.
//

#include "File.h"


bool hdf5::FileBase::is_hdf5(const std::string &fname) {
    return H5Fis_hdf5(fname.c_str());
}

void hdf5::FileBase::require_is_hdf5(const std::string &fname) {
    REQUIRE_TRUE(is_hdf5(fname), "Specified file is not HDF5 format");
}