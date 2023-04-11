//
// Created by anderson on 27/06/2022.
//

#include "File.h"

bool hdf5::FileBase::is_hdf5(const str_t& fname) {
    if (fname.empty()) return false;
    return H5Fis_hdf5(fname.c_str()) > 0;
}

void hdf5::FileBase::require_is_hdf5(const str_t& fname) {
    REQUIRE_TRUE(is_hdf5(fname), "Specified file is not HDF5 format");
}

hid_t hdf5::FileReader::get_id(const str_t& fname) {
    require_is_hdf5(fname);
    auto plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    REQUIRE_TRUE(H5Fis_hdf5(fname.c_str()), "Specified file is not HDF5 format");
    auto id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist);
    H5Pclose(plist);
    return id;
}

hid_t hdf5::FileWriter::get_id(const str_t& fname) {
    auto plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    auto id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    REQUIRE_NE(id, 0, "HDF5 file could not be opened for writing. It may be locked by another program");
    H5Pclose(plist);
    return id;
}
