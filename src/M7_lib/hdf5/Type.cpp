//
// Created by anderson on 27/06/2022.
//

#include "Type.h"

hsize_t hdf5::type_size(hid_t h5type) {
    return H5Tget_size(h5type);
}


hsize_t hdf5::StringType::size_max(const std::vector<std::string>& vec) {
    if (vec.empty()) return 0ul;
    return std::max_element(vec.cbegin(), vec.cend(),
                            [](const std::string& s1, const std::string& s2){return s1.size()>s2.size();})->size();
}

hdf5::StringType::StringType(hsize_t size, int) : m_handle(H5Tcopy(H5T_C_S1)), m_nchar(size+1) {
    auto status = H5Tset_size(m_handle, m_nchar);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type resizing failed");
    DEBUG_ASSERT_EQ(H5Tget_size(m_handle), m_nchar, "string length at odds with type length");
}

hdf5::StringType::StringType(hid_t handle) : m_handle(handle), m_nchar(H5Tget_size(m_handle)){
    DEBUG_ASSERT_TRUE(m_nchar, "number of chars in string type should be non-zero");
}

hdf5::StringType::StringType(const std::string& str) : StringType(str.size(), 0) {}

hdf5::StringType::StringType(const std::vector<std::string>& str_vec) : StringType(size_max(str_vec), 0) {}

hdf5::StringType::~StringType() {
    auto status = H5Tclose(m_handle);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type release failed");
}

hdf5::StringType::operator hid_t() const {
    return m_handle;
}
