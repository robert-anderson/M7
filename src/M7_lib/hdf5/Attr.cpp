//
// Created by anderson on 27/06/2022.
//

#include "Attr.h"


hdf5::AttrReader::AttrReader(hid_t parent_handle, const str_t& name) :
        m_handle(H5Aopen(parent_handle, name.c_str(), H5P_DEFAULT)),
        m_space(H5Aget_space(m_handle)), m_type(H5Aget_type(m_handle)),
        m_nelement(m_space.m_nelement){}

hdf5::AttrReader::~AttrReader() {
    H5Aclose(m_handle);
}

void hdf5::AttrReader::read_bytes(char* dst) const {
    auto status = H5Aread(m_handle, m_type, dst);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute read failed");
}

void hdf5::AttrReader::read(str_t* dst, size_t n) const {
    REQUIRE_EQ(n, m_space.m_nelement, "number of elements read must be the number stored");
    std::vector<char> tmp(m_type.m_size*n);
    auto status = H5Aread(m_handle, m_type, tmp.data());
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute read failed");
    for (uint_t i=0; i<n; ++i) {
        (dst++)->insert(0, tmp.data()+i*m_type.m_size, m_type.m_size);
    }
}

hdf5::AttrWriter::AttrWriter(hid_t parent_handle, const str_t& name, const std::vector<hsize_t>& shape,
                             hid_t h5type) :
        m_space(shape), m_h5type(h5type),
        m_handle(H5Acreate(parent_handle, name.c_str(), m_h5type, m_space.m_handle, H5P_DEFAULT, H5P_DEFAULT)){}

hdf5::AttrWriter::~AttrWriter() {
    H5Aclose(m_handle);
}

void hdf5::AttrWriter::write_bytes(const char* src) const {
    auto status = H5Awrite(m_handle, m_h5type, src);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute write failed");
}