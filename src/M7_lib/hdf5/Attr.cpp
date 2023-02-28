//
// Created by anderson on 27/06/2022.
//

#include "Attr.h"


hdf5::Attr::Attr(v_t<buf_t> buf, hdf5::dataset::ItemFormat format, str_t name, char) :
        m_buf(std::move(buf)), m_format(std::move(format)), m_name(std::move(name)){
    REQUIRE_EQ(m_buf.size(), m_format.m_size, "buffer size inconsistent with format");
}

hdf5::Attr::Attr(hid_t parent_id, str_t name) : Attr(load(parent_id, name)){}

bool hdf5::Attr::operator==(const hdf5::Attr& other) const {
    if (m_buf != other.m_buf) return false;
    if (m_format!=other.m_format) return false;
    if (m_name!=other.m_name) return false;
    return true;
}

bool hdf5::Attr::operator!=(const hdf5::Attr& other) const {
    return !(*this==other);
}

void hdf5::Attr::save(hid_t parent_id) const {
    auto dataspace_id = H5Screate_simple(m_format.m_h5_shape.size(), m_format.m_h5_shape.data(), nullptr);
    auto attr_id = H5Acreate(parent_id, m_name.c_str(), m_format.m_type.m_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    auto status = H5Awrite(attr_id, m_format.m_type.m_id, m_buf.data());
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute write failed");
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);
}

hdf5::Attr hdf5::Attr::load(hid_t parent_id, const str_t& name) {
    if (!H5Aexists(parent_id, name.c_str())) return {{}, {}, name, 0};
    auto attr_id = H5Aopen(parent_id, name.c_str(), H5P_DEFAULT);
    auto dataspace_id = H5Aget_space(attr_id);
    auto ndim = H5Sget_simple_extent_ndims(dataspace_id);
    v_t<hsize_t> shape(ndim);
    H5Sget_simple_extent_dims(dataspace_id, shape.data(), nullptr);
    Type type(H5Aget_type(attr_id));
    dataset::ItemFormat format(type, convert::vector<uint_t>(shape), {}, false);
    v_t<buf_t> buf(format.m_size);
    auto status = H5Aread(attr_id, type.m_id, buf.data());
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute write failed");
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);
    return {buf, format, name, 0};
}
