//
// Created by rja on 13/12/2020.
//

#include "HDF5Wrapper.h"

hdf5::Group hdf5::File::subgroup(std::string name) {
    return Group(m_handle, name, m_writemode);
}

hdf5::AttributeWriterBase::AttributeWriterBase(hid_t parent_handle, std::string name, const defs::inds &shape,
                                               hid_t h5type) :
        m_parent_handle(parent_handle), m_h5type(h5type), m_shape(shape),
        m_nelement(nd_utils::nelement(shape)) {
    auto shape_tmp = convert_dims(shape);
    m_memspace_handle = H5Screate_simple(shape.size(), shape_tmp.data(), nullptr);
    m_handle = H5Acreate(m_parent_handle, name.c_str(), h5type, m_memspace_handle, H5P_DEFAULT, H5P_DEFAULT);
}

hdf5::AttributeWriterBase::~AttributeWriterBase() {
    H5Sclose(m_memspace_handle);
    H5Aclose(m_handle);
}

void hdf5::AttributeWriterBase::write_bytes(const char *src) {
    auto status = H5Awrite(m_handle, m_h5type, src);
    MPI_ASSERT(!status, "HDF5 attribute write failed");
}

void hdf5::AttributeWriterBase::write(hid_t parent, std::string name, const std::string &src) {
    auto type = H5Tcopy(H5T_C_S1);
    auto status = H5Tset_size(type, src.size() + 1); // include space for null terminator
    MPI_ASSERT(!status, "HDF5 string type resizing failed");
    AttributeWriterBase(parent, name, {1}, type).write_bytes((const char *) src.c_str());
    status = H5Tclose(type);
    MPI_ASSERT(!status, "HDF5 string type release failed");
}

void hdf5::AttributeWriterBase::write(hid_t parent, std::string name, const std::vector<std::string> &src) {
    // find longest string in vector
    size_t max_size=0ul;
    for (const auto& str: src) max_size = (str.size()>max_size) ? str.size() : max_size;
    max_size++; // include space for null terminator
    std::string tmp(max_size*src.size(), 0);
    auto tmp_it = tmp.begin();
    for (const auto& str: src) {
        std::copy(str.cbegin(), str.cend(), tmp_it);
        tmp_it+=max_size;
    }
    // now build the type
    auto type = H5Tcopy(H5T_C_S1);
    auto status = H5Tset_size(type, max_size);
    MPI_ASSERT(!status, "HDF5 string type resizing failed");
    AttributeWriterBase(parent, name, {src.size()}, type).write_bytes((const char *) tmp.c_str());
    status = H5Tclose(type);
    MPI_ASSERT(!status, "HDF5 string type release failed");
}
