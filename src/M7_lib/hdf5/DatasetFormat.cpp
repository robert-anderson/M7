//
// Created by rja on 06/12/22.
//

#include "DatasetFormat.h"
#include "M7_lib/util/Vector.h"

hdf5::dataset::ItemFormat::ItemFormat(hdf5::Type type, uintv_t shape, strv_t dim_names, bool add_complex_dim) :
        m_type(type), m_shape(std::move(shape)),
        m_h5_shape(convert::vector<hsize_t>(add_complex_dim ? vector::appended(m_shape, 2ul) : m_shape)),
        m_size(nd::nelement(m_h5_shape) * m_type.m_size),
        m_dim_names(add_complex_dim && !dim_names.empty() ?
                    vector::appended(dim_names, "real/imag") : dim_names) {
    REQUIRE_TRUE(m_dim_names.empty() || m_dim_names.size() == m_h5_shape.size(),
                 "incorrect number of dimension names");
}

bool hdf5::dataset::ItemFormat::operator==(const hdf5::dataset::ItemFormat& other) const {
    return m_type==other.m_type && m_h5_shape==other.m_h5_shape && m_dim_names==other.m_dim_names;
}

bool hdf5::dataset::ItemFormat::operator!=(const hdf5::dataset::ItemFormat& other) const {
    return !(*this==other);
}

uint_t hdf5::dataset::ListFormat::ndim() const {
    return m_h5_shape.size();
}

hdf5::dataset::ListFormat::ListFormat(hdf5::dataset::ItemFormat item_format, uint_t nitem, str_t leading_dim_name) :
        m_item(item_format), m_nitem(nitem), m_size(m_nitem*m_item.m_size),
        m_h5_shape(convert::vector<hsize_t>(vector::prepended(m_item.m_h5_shape, m_nitem))),
        m_dim_names(vector::prepended(m_item.m_dim_names, leading_dim_name)) {}


hdf5::dataset::DistListFormat::DistListFormat(
        hdf5::dataset::ItemFormat item_format, uint_t nitem_local,
        uint_t nitem, uint_t nitem_displ, str_t leading_dim_name) :
        m_local(item_format, nitem_local, leading_dim_name), m_nitem(nitem), m_nitem_displ(nitem_displ),
        m_h5_shape(vector::prepended(m_local.m_item.m_h5_shape, m_nitem)){}

hdf5::dataset::PartDistListFormat::PartDistListFormat(
        hdf5::dataset::ItemFormat item_format, uint_t nitem, str_t leading_dim_name) :
        DistListFormat(item_format, nitem, mpi::all_sum(nitem),
            mpi::counts_to_displs_consec(mpi::all_gathered(nitem))[mpi::irank()], leading_dim_name){}

hdf5::dataset::FullDistListFormat::FullDistListFormat(
        hdf5::dataset::ItemFormat item_format, uint_t nitem, str_t leading_dim_name) :
        DistListFormat(item_format, nitem, mpi::all_max(nitem), 0ul, leading_dim_name){
    if (m_local.m_nitem) REQUIRE_EQ(m_nitem, m_local.m_nitem,
                            "full distributed list requires that participating ranks read all items");
}
