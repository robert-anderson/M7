//
// Created by anderson on 27/06/2022.
//

#include "Type.h"

hsize_t hdf5::type_size(hid_t h5type) {
    return H5Tget_size(h5type);
}


hsize_t hdf5::Type::size_max(const std::vector<std::string>* vec) {
    if (!vec || vec->empty()) return 0ul;
    return std::max_element(vec->cbegin(), vec->cend(),
                            [](const std::string& s1, const std::string& s2){return s1.size()>s2.size();})->size();
}

hdf5::Type::Type(hsize_t size, char) : m_handle(H5Tcopy(H5T_C_S1)), m_size(size), m_immutable(false){
    auto status = H5Tset_size(m_handle, size+1);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type resizing failed");
    DEBUG_ASSERT_EQ(H5Tget_size(m_handle), m_size+1, "string length at odds with type length");
}

hdf5::Type::Type(const std::string* str) : Type(str->size(), 0) {}

hdf5::Type::Type(const std::vector<std::string>* str_vec) : Type(size_max(str_vec), 0) {}

hdf5::Type::~Type() {
    if (m_immutable) return;
    auto status = H5Tclose(m_handle);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type release failed");
}

hdf5::Type::operator hid_t() const {
    return m_handle;
}
