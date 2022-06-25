//
// Created by Robert J. Anderson on 13/12/2020.
//

#include "HDF5Wrapper.h"



hsize_t hdf5::type_size(hid_t h5type) {
    return H5Tget_size(h5type);
}


void hdf5::FileBase::check_is_hdf5(const std::string &name) {
    REQUIRE_TRUE(H5Fis_hdf5(name.c_str()), "Specified file is not HDF5 format");
}

bool hdf5::NodeReader::child_exists(const std::string &name) const {
    for (uint_t i=0ul; i<nchild(); ++i) if (name==child_name(i)) return true;
    return false;
}

uint_t hdf5::NodeReader::first_existing_child(const std::vector<std::string> &names) const {
    for (uint_t i=0ul; i<names.size(); ++i) if (child_exists(names[i])) return i;
    return ~0ul;
}

uint_t hdf5::NodeReader::nchild() const {
    hsize_t num;
    auto status = H5Gget_num_objs(m_handle, &num);
    REQUIRE_TRUE(!status, "could not get number of objects within HDF5 group");
    return status == 0;
}

std::string hdf5::NodeReader::child_name(uint_t ichild) const {
    uint_t size = H5Lget_name_by_idx(m_handle, ".", H5_INDEX_NAME,
                                     H5_ITER_INC, ichild, nullptr, 0, H5P_DEFAULT);
    std::string name(size, 0);
    auto name_ptr = const_cast<char*>(name.c_str());
    uint_t size_chk = H5Lget_name_by_idx(m_handle, ".", H5_INDEX_NAME,
                                         H5_ITER_INC, ichild, name_ptr, size+1, H5P_DEFAULT);
    REQUIRE_EQ(size, size_chk, "inconsistent number of chars in name");
    return name;
}

int hdf5::NodeReader::child_type(uint_t i) const {
    return H5Gget_objtype_by_idx(m_handle, i);
}

std::vector<std::string> hdf5::NodeReader::child_names(int type) const {
    std::vector<std::string> names;
    auto n = nchild();
    names.reserve(n);
    for (uint_t i=0ul; i<n; ++i) {
        if (type<0 || child_type(i)==type) names.push_back(child_name(i));
    }
    return names;
}

uint_t hdf5::NodeReader::get_dataset_ndim(std::string name) const {
    auto status = H5Gget_objinfo(m_handle, name.c_str(), 0, nullptr);
    REQUIRE_TRUE(!status, "Dataset \"" + name + "\" does not exist");
    auto dataset = H5Dopen1(m_handle, name.c_str());
    auto dataspace = H5Dget_space(dataset);
    auto rank = H5Sget_simple_extent_ndims(dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return rank;
}

uintv_t hdf5::NodeReader::get_dataset_shape(std::string name) const {
    auto ndim = get_dataset_ndim(name);
    auto dataset = H5Dopen1(m_handle, name.c_str());
    REQUIRE_GT_ALL(dataset, 0, log::format("no such dataset \"{}\"", name));
    auto dataspace = H5Dget_space(dataset);
    std::vector<hsize_t> dims(ndim, 0ul);
    H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    uintv_t out;
    out.reserve(dims.size());
    for (const auto &i: dims) out.push_back(i);
    return out;
}


void hdf5::NdDistListWriter::write_h5item_bytes(const uint_t &iitem, const void *data) {
    DEBUG_ASSERT_EQ(bool(data), iitem < m_nitem_local,
               "data is null and items remain, or this is a runoff write op and data is not null");
    select_hyperslab(iitem);
    log::debug_("writing data...");

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
    log::debug_("data written");
}

hdf5::NdDistListWriter::NdDistListWriter(hid_t parent_handle, std::string name, const uintv_t &item_dims,
                                         const uint_t &nitem, hid_t h5type, const std::vector<std::string> &dim_labels)
        :NdDistListBase(parent_handle, name, item_dims, nitem, true, h5type) {
    if (!dim_labels.empty()) {
        DEBUG_ASSERT_EQ(dim_labels.size(), item_dims.size(),
                        "Number of dim labels does not match number of dims");
        for (uint_t idim = 0ul; idim < item_dims.size(); ++idim)
            H5DSset_label(m_dataset_handle, idim, dim_labels[idim].c_str());
    }
}

std::vector<hsize_t> hdf5::NdDistListBase::get_list_dims_local() {
    std::vector<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_local);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

std::vector<hsize_t> hdf5::NdDistListBase::get_list_dims_global() {
    std::vector<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_global);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

hsize_t hdf5::NdDistListBase::get_item_offset() {
    std::vector<hsize_t> tmp(mpi::nrank());
    mpi::all_gather(m_nitem_local, tmp);
    hsize_t out = 0ul;
    for (uint_t irank = 0ul; irank < mpi::irank(); ++irank) out += tmp[irank];
    return out;
}

hdf5::NdDistListBase::NdDistListBase(hid_t parent_handle, std::string name, const uintv_t &item_dims, const uint_t &nitem,
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
    std::vector<hsize_t> zeros(m_ndim_list, 0ul);
    m_none_memspace_handle = H5Screate_simple(m_ndim_list, zeros.data(), nullptr);

    log::debug_("Opened HDF5 NdList with {} local items", m_nitem_local);
}

void hdf5::NdDistListBase::select_hyperslab(const uint_t &iitem) {
    if (iitem < m_nitem_local) {
        m_hyperslab_offsets[0] = m_item_offset + iitem;
        log::debug_("selecting hyperslab with offsets: {}", convert::to_string(m_hyperslab_offsets));
        auto status = H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                          nullptr, m_hyperslab_counts.data(), nullptr);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 hyperslab selection failed");
        log::debug_("hyperslab selected");
    } else {
        /*
         * the item index exceeds local bounds, so make a null selection instead
         * of a hyperslab selection
         */
        log::debug_("making null selection");
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
