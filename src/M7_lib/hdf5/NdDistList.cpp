//
// Created by anderson on 27/06/2022.
//

#include "NdDistList.h"


v_t<hsize_t> hdf5::DistBase::get_list_dims_local() {
    v_t<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_local);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

v_t<hsize_t> hdf5::DistBase::get_list_dims_global() {
    v_t<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_global);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

hsize_t hdf5::DistBase::get_item_offset() {
    v_t<hsize_t> tmp(mpi::nrank());
    mpi::all_gather(m_nitem_local, tmp);
    hsize_t out = 0ul;
    for (uint_t irank = 0ul; irank < mpi::irank(); ++irank) out += tmp[irank];
    return out;
}

hdf5::DistBase::DistBase(hid_t parent_handle, str_t name, const uintv_t &item_dims, const uint_t &nitem,
                                     bool writemode, hid_t h5type) :
        m_parent_handle(parent_handle),
        m_item_dims(convert::vector<hsize_t>(item_dims)),
        m_ndim_item(item_dims.size()),
        m_ndim_list(item_dims.size() + 1),
        m_nitem_local(nitem),
        m_nitem_global(mpi::all_sum(m_nitem_local)),
        m_nitem_local_max(mpi::all_max(m_nitem_local)),
        m_list_dims_local(get_list_dims_local()),
        m_list_dims_global(get_list_dims_global()),
        m_item_offset(get_item_offset()),
        m_hyperslab_counts(m_ndim_list, 0ul),
        m_hyperslab_offsets(m_ndim_list, 0ul),
        m_h5type(h5type) {
    REQUIRE_TRUE(H5Tget_size(m_h5type), "Invalid HDF5 type specified");
    m_filespace_handle = H5Screate_simple(m_ndim_list, m_list_dims_global.data(), nullptr);

    /*
     * Create the dataset with default properties and close filespace.
     */
    if (writemode)
        m_dataset_handle = H5Dcreate(m_parent_handle, name.c_str(), m_h5type, m_filespace_handle,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    else
        m_dataset_handle = H5Dopen1(m_parent_handle, name.c_str());

    H5Sclose(m_filespace_handle);

    m_filespace_handle = H5Dget_space(m_dataset_handle);

    m_hyperslab_counts = m_list_dims_local;
    // select one item at a time
    m_hyperslab_counts[0] = 1;
    m_memspace_handle = H5Screate_simple(m_ndim_list, m_hyperslab_counts.data(), nullptr);
    v_t<hsize_t> zeros(m_ndim_list, 0ul);
    m_none_memspace_handle = H5Screate_simple(m_ndim_list, zeros.data(), nullptr);

    logging::debug_("Opened HDF5 NdList with {} local items", m_nitem_local);
}

void hdf5::DistBase::select_hyperslab(const uint_t &iitem) {
    if (iitem < m_nitem_local) {
        m_hyperslab_offsets[0] = m_item_offset + iitem;
        logging::debug_("selecting hyperslab with offsets: {}", convert::to_string(m_hyperslab_offsets));
        auto status = H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                          nullptr, m_hyperslab_counts.data(), nullptr);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 hyperslab selection failed");
        logging::debug_("hyperslab selected");
    } else {
        /*
         * the item index exceeds local bounds, so make a null selection instead
         * of a hyperslab selection
         */
        logging::debug_("making null selection");
        DEBUG_ASSERT_LT(iitem, m_nitem_local_max, "Item index exceeds global maximum");
        auto status = H5Sselect_none(m_filespace_handle);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 null selection failed");
    }
}

void hdf5::DistBase::select_hyperslab(uint_t iitem_begin, uint_t nitem) {
    DEBUG_ASSERT_LE(iitem_begin+nitem, m_nitem_local, "last item OOB");
    if (nitem) {
        m_hyperslab_offsets[0] = m_item_offset + iitem_begin;
        m_hyperslab_counts[0] = nitem;
        logging::debug_("selecting hyperslab with offsets: {}", convert::to_string(m_hyperslab_offsets));
        auto status = H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                          nullptr, m_hyperslab_counts.data(), nullptr);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 hyperslab selection failed");
        logging::debug_("hyperslab selected");
    } else {
        /*
         * nothing to write, so make a null selection instead of a hyperslab selection
         */
        logging::debug_("making null selection");
        DEBUG_ASSERT_LT(iitem, m_nitem_local_max, "Item index exceeds global maximum");
        auto status = H5Sselect_none(m_filespace_handle);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 null selection failed");
    }
}


hdf5::DistBase::~DistBase() {
    H5Sclose(m_filespace_handle);
    H5Dclose(m_dataset_handle);
    H5Sclose(m_memspace_handle);
}


















void hdf5::NdDistListWriter::write_h5item_bytes(const uint_t &iitem, const void *data) {
    DEBUG_ASSERT_EQ(bool(data), iitem < m_nitem_local,
                    "data is null and items remain, or this is a runoff write op and data is not null");
    select_hyperslab(iitem);
    logging::debug_("writing data...");

    if (data) {
        auto status = H5Dwrite(m_dataset_handle, m_h5type, m_memspace_handle,
                               m_filespace_handle, m_coll_plist, data);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 write failed");
    }
    else {
        auto status = H5Dwrite(m_dataset_handle, m_h5type, m_none_memspace_handle,
                               m_filespace_handle, m_coll_plist, data);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 write failed");
    }
    logging::debug_("data written");
}


void hdf5::NdDistListWriter::write(const void *data, uint_t iitem_begin, uint_t nitem) {
    DEBUG_ASSERT_EQ(bool(data), bool(nitem),
                    "data is null and items remain, or this is a runoff write op and data is not null");
    select_hyperslab(iitem_begin, nitem);
    logging::debug_("writing data...");

    if (data) {
        auto status = H5Dwrite(m_dataset_handle, m_h5type, m_memspace_handle,
                               m_filespace_handle, m_coll_plist, data);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 write failed");
    }
    else {
        auto status = H5Dwrite(m_dataset_handle, m_h5type, m_none_memspace_handle,
                               m_filespace_handle, m_coll_plist, data);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 write failed");
    }
    logging::debug_("data written");
}

hdf5::NdDistListWriter::NdDistListWriter(hid_t parent_handle, str_t name, const uintv_t &item_dims,
                                         const uint_t &nitem, hid_t h5type, const strv_t &dim_labels)
        :NdDistListBase(parent_handle, name, item_dims, nitem, true, h5type) {
    if (!dim_labels.empty()) {
        DEBUG_ASSERT_EQ(dim_labels.size(), item_dims.size(),
                        "Number of dim labels does not match number of dims");
        for (uint_t idim = 0ul; idim < item_dims.size(); ++idim)
            H5DSset_label(m_dataset_handle, idim, dim_labels[idim].c_str());
    }
}

v_t<hsize_t> hdf5::NdDistListBase::get_list_dims_local() {
    v_t<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_local);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

v_t<hsize_t> hdf5::NdDistListBase::get_list_dims_global() {
    v_t<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_global);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

hsize_t hdf5::NdDistListBase::get_item_offset() {
    v_t<hsize_t> tmp(mpi::nrank());
    mpi::all_gather(m_nitem_local, tmp);
    hsize_t out = 0ul;
    for (uint_t irank = 0ul; irank < mpi::irank(); ++irank) out += tmp[irank];
    return out;
}

hdf5::NdDistListBase::NdDistListBase(hid_t parent_handle, str_t name, const uintv_t &item_dims, const uint_t &nitem,
                                     bool writemode, hid_t h5type) :
        m_parent_handle(parent_handle),
        m_item_dims(convert::vector<hsize_t>(item_dims)),
        m_ndim_item(item_dims.size()),
        m_ndim_list(item_dims.size() + 1),
        m_nitem_local(nitem),
        m_nitem_global(mpi::all_sum(m_nitem_local)),
        m_nitem_local_max(mpi::all_max(m_nitem_local)),
        m_list_dims_local(get_list_dims_local()),
        m_list_dims_global(get_list_dims_global()),
        m_item_offset(get_item_offset()),
        m_hyperslab_counts(m_ndim_list, 0ul),
        m_hyperslab_offsets(m_ndim_list, 0ul),
        m_h5type(h5type) {
    REQUIRE_TRUE(H5Tget_size(m_h5type), "Invalid HDF5 type specified");
    m_filespace_handle = H5Screate_simple(m_ndim_list, m_list_dims_global.data(), nullptr);

    /*
     * Create the dataset with default properties and close filespace.
     */
    if (writemode)
        m_dataset_handle = H5Dcreate(m_parent_handle, name.c_str(), m_h5type, m_filespace_handle,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    else
        m_dataset_handle = H5Dopen1(m_parent_handle, name.c_str());

    H5Sclose(m_filespace_handle);

    m_filespace_handle = H5Dget_space(m_dataset_handle);

    m_hyperslab_counts = m_list_dims_local;
    // select one item at a time
    m_hyperslab_counts[0] = 1;
    m_memspace_handle = H5Screate_simple(m_ndim_list, m_hyperslab_counts.data(), nullptr);
    v_t<hsize_t> zeros(m_ndim_list, 0ul);
    m_none_memspace_handle = H5Screate_simple(m_ndim_list, zeros.data(), nullptr);

    logging::debug_("Opened HDF5 NdList with {} local items", m_nitem_local);
}

void hdf5::NdDistListBase::select_hyperslab(const uint_t &iitem) {
    if (iitem < m_nitem_local) {
        m_hyperslab_offsets[0] = m_item_offset + iitem;
        logging::debug_("selecting hyperslab with offsets: {}", convert::to_string(m_hyperslab_offsets));
        auto status = H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                          nullptr, m_hyperslab_counts.data(), nullptr);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 hyperslab selection failed");
        logging::debug_("hyperslab selected");
    } else {
        /*
         * the item index exceeds local bounds, so make a null selection instead
         * of a hyperslab selection
         */
        logging::debug_("making null selection");
        DEBUG_ASSERT_LT(iitem, m_nitem_local_max, "Item index exceeds global maximum");
        auto status = H5Sselect_none(m_filespace_handle);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 null selection failed");
    }
}

void hdf5::NdDistListBase::select_hyperslab(uint_t iitem_begin, uint_t nitem) {
    DEBUG_ASSERT_LE(iitem_begin+nitem, m_nitem_local, "last item OOB");
    if (nitem) {
        m_hyperslab_offsets[0] = m_item_offset + iitem_begin;
        m_hyperslab_counts[0] = nitem;
        logging::debug_("selecting hyperslab with offsets: {}", convert::to_string(m_hyperslab_offsets));
        auto status = H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                          nullptr, m_hyperslab_counts.data(), nullptr);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 hyperslab selection failed");
        logging::debug_("hyperslab selected");
    } else {
        /*
         * nothing to write, so make a null selection instead of a hyperslab selection
         */
        logging::debug_("making null selection");
        DEBUG_ASSERT_LT(iitem, m_nitem_local_max, "Item index exceeds global maximum");
        auto status = H5Sselect_none(m_filespace_handle);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 null selection failed");
    }
}


hdf5::NdDistListBase::~NdDistListBase() {
    H5Sclose(m_filespace_handle);
    H5Dclose(m_dataset_handle);
    H5Sclose(m_memspace_handle);
}
