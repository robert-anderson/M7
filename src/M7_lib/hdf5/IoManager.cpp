//
// Created by rja on 06/12/22.
//

#include "IoManager.h"
#include "M7_lib/util/Vector.h"

hdf5::dataset::ItemFormat::ItemFormat(hdf5::Type type, uintv_t shape, strv_t dim_names, bool add_complex_dim) :
        m_type(type), m_shape(std::move(shape)),
        m_h5_shape(convert::vector<hsize_t>(add_complex_dim ? vector::appended(m_shape, 2ul) : m_shape)),
        m_size(nd::nelement(m_h5_shape) * m_type.m_size),
        m_dim_names(add_complex_dim && !dim_names.empty() ?
                    vector::appended(dim_names, "real/imag") : dim_names) {
    REQUIRE_TRUE(m_dim_names.empty() || m_dim_names.size() == m_h5_shape.size(),
                 "incorrect number of dimension names");
    REQUIRE_TRUE(m_size, "format is empty");
}

bool hdf5::dataset::ItemFormat::operator==(const hdf5::dataset::ItemFormat& other) const {
    return m_type==other.m_type && m_h5_shape==other.m_h5_shape && m_dim_names==other.m_dim_names;
}

bool hdf5::dataset::ItemFormat::operator!=(const hdf5::dataset::ItemFormat& other) const {
    return !(*this==other);
}

hdf5::dataset::ListFormat::ListFormat(hdf5::dataset::ItemFormat item_format, uint_t nitem) :
        m_item(item_format), m_nitem(nitem), m_size(m_nitem*m_item.m_size),
        m_h5_shape(convert::vector<hsize_t>(vector::prepended(m_item.m_h5_shape, m_nitem))),
        m_dim_names(vector::prepended(m_item.m_dim_names, "item")) {}

hdf5::dataset::DistListFormat::DistListFormat(
        hdf5::dataset::ItemFormat item_format, uint_t nitem_local, uint_t nitem, uint_t nitem_displ) :
        m_local(item_format, nitem_local), m_nitem(nitem), m_nitem_displ(nitem_displ),
        m_h5_shape(vector::prepended(m_local.m_item.m_h5_shape, m_nitem)){}
