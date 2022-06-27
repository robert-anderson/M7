//
// Created by anderson on 27/06/2022.
//

#include "Dataspace.h"

std::vector<hsize_t> hdf5::DataSpace::make_shape() const {
    auto ndim = H5Sget_simple_extent_dims(m_handle, nullptr, nullptr);
    std::vector<hsize_t> shape(ndim, 0);
    H5Sget_simple_extent_dims(m_handle, shape.data(), nullptr);
    return shape;
}

hdf5::DataSpace::DataSpace(hid_t handle) : m_handle(handle), m_shape(make_shape()), m_nelement(nd::nelement(m_shape)){}

hdf5::DataSpace::DataSpace(const std::vector<hsize_t>& shape, bool select_none) :
        DataSpace(H5Screate_simple(shape.size(), shape.data(), nullptr)){
    REQUIRE_EQ(shape, m_shape, "given shape and shape reported by HDF5 do not agree");
    if (select_none) this->select_none();
}

hdf5::DataSpace::~DataSpace() {
    H5Sclose(m_handle);
}

hdf5::DataSpace::operator hid_t() const {
    return m_handle;
}

void hdf5::DataSpace::select_none() const {
    H5Sselect_none(m_handle);
}
