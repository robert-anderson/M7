//
// Created by anderson on 27/06/2022.
//

#include "Group.h"

hdf5::GroupWriter::GroupWriter(const hdf5::NodeWriter &node, str_t name) :
        NodeWriter(H5Gcreate(node, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {}

hdf5::GroupWriter::~GroupWriter() {
    auto status = H5Gclose(m_handle);
    REQUIRE_TRUE(!status, "HDF5 Error on closing group");
}

hdf5::GroupReader::GroupReader(const hdf5::NodeReader &node, str_t name) :
        NodeReader(H5Gopen2(node, name.c_str(), H5P_DEFAULT)) {}

hdf5::GroupReader::~GroupReader() {
    auto status = H5Gclose(m_handle);
    REQUIRE_TRUE(!status, "HDF5 Error on closing group");
}
