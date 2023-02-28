//
// Created by anderson on 27/06/2022.
//

#include "Type.h"

hsize_t hdf5::Type::size_max(const strv_t* vec) {
    if (!vec || vec->empty()) return 0ul;
    return std::max_element(vec->cbegin(), vec->cend(),
                            [](const str_t& s1, const str_t& s2){return s1.size()>s2.size();})->size();
}

hdf5::Type::Type(hsize_t size, tag::Int<true> /*is_string*/) :
        m_id(H5Tcopy(H5T_C_S1)), m_size(size), m_immutable(false){
    auto status = H5Tset_size(m_id, size);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type resizing failed");
    DEBUG_ASSERT_EQ(H5Tget_size(m_id), m_size, "string length at odds with type length");
}

hdf5::Type::Type(const str_t* str) : Type(str->size(), tag::Int<true>()) {}

hdf5::Type::Type(const strv_t* strv) : Type(size_max(strv), tag::Int<true>()) {}

hdf5::Type::Type(const hdf5::Type& other) :
        m_id(other.m_immutable ? other.m_id : H5Tcopy(other.m_id)),
        m_size(other.m_size), m_immutable(other.m_immutable){}

hdf5::Type::~Type() {
    if (m_immutable) return;
    auto status = H5Tclose(m_id);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type release failed");
}