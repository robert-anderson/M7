//
// Created by rja on 13/12/2020.
//

#include "HDF5Wrapper.h"

hdf5::Group hdf5::File::subgroup(std::string name) {
    return Group(m_handle, name, m_writemode);
}

hdf5::AttributeWriterBase::AttributeWriterBase(hid_t parent_handle, std::string name, const defs::inds &shape,
                                               hid_t h5type) :
        m_parent_handle(parent_handle), m_h5type(h5type), m_shape(shape),
        m_nelement(nd_utils::nelement(shape)) {
    auto shape_tmp = convert_dims(shape);
    m_memspace_handle = H5Screate_simple(shape.size(), shape_tmp.data(), nullptr);
    m_handle = H5Acreate(m_parent_handle, name.c_str(), h5type, m_memspace_handle, H5P_DEFAULT, H5P_DEFAULT);
}

hdf5::AttributeWriterBase::~AttributeWriterBase() {
    H5Sclose(m_memspace_handle);
    H5Aclose(m_handle);
}

void hdf5::AttributeWriterBase::write_bytes(const char *src) {
    auto status = H5Awrite(m_handle, m_h5type, src);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute write failed");
}

void hdf5::AttributeWriterBase::write(hid_t parent, std::string name, const std::string &src) {
    auto type = H5Tcopy(H5T_C_S1);
    auto status = H5Tset_size(type, src.size() + 1); // include space for null terminator
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type resizing failed");
    AttributeWriterBase(parent, name, {1}, type).write_bytes((const char *) src.c_str());
    status = H5Tclose(type);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type release failed");
}

void hdf5::AttributeWriterBase::write(hid_t parent, std::string name, const std::vector<std::string> &src) {
    // find longest string in vector
    size_t max_size = 0ul;
    for (const auto &str: src) max_size = (str.size() > max_size) ? str.size() : max_size;
    max_size++; // include space for null terminator
    std::string tmp(max_size * src.size(), 0);
    auto tmp_it = tmp.begin();
    for (const auto &str: src) {
        std::copy(str.cbegin(), str.cend(), tmp_it);
        tmp_it += max_size;
    }
    // now build the type
    auto type = H5Tcopy(H5T_C_S1);
    auto status = H5Tset_size(type, max_size);
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type resizing failed");
    AttributeWriterBase(parent, name, {src.size()}, type).write_bytes((const char *) tmp.c_str());
    status = H5Tclose(type);
    DEBUG_ASSERT_FALSE(status, "HDF5 string type release failed");
}

void hdf5::NdListWriter::write_h5item_bytes(const size_t &iitem, const void *data) {
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

hdf5::NdListWriter::NdListWriter(hid_t parent_handle, std::string name, const defs::inds &item_dims,
                                 const size_t &nitem, hid_t h5type, const std::vector<std::string> &dim_labels)
        :
        NdListBase(parent_handle, name, item_dims, nitem, true, h5type) {
    if (!dim_labels.empty()) {
        DEBUG_ASSERT_EQ(dim_labels.size(), item_dims.size(),
                        "Number of dim labels does not match number of dims");
        for (size_t idim = 0ul; idim < item_dims.size(); ++idim)
            H5DSset_label(m_dataset_handle, idim, dim_labels[idim].c_str());
    }
}

std::vector<hsize_t> hdf5::NdListBase::get_list_dims_local() {
    std::vector<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_local);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

std::vector<hsize_t> hdf5::NdListBase::get_list_dims_global() {
    std::vector<hsize_t> out;
    out.reserve(m_ndim_list);
    out.push_back(m_nitem_global);
    out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
    return out;
}

hsize_t hdf5::NdListBase::get_item_offset() {
    std::vector<hsize_t> tmp(mpi::nrank());
    mpi::all_gather(m_nitem_local, tmp);
    hsize_t out = 0ul;
    for (size_t irank = 0ul; irank < mpi::irank(); ++irank) out += tmp[irank];
    return out;
}

hdf5::NdListBase::NdListBase(hid_t parent_handle, std::string name, const defs::inds &item_dims, const size_t &nitem,
                             bool writemode, hid_t h5type) :
        m_parent_handle(parent_handle),
        m_item_dims(convert_dims(item_dims)),
        m_ndim_item(item_dims.size()),
        m_ndim_list(item_dims.size() + 1),
        m_nitem_local(nitem),
        m_nitem_global(mpi::all_sum(m_nitem_local)),
        m_nitem_global_max(mpi::all_max(m_nitem_local)),
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

void hdf5::NdListBase::select_hyperslab(const size_t &iitem) {
    if (iitem < m_nitem_local) {
        m_hyperslab_offsets[0] = m_item_offset + iitem;
        log::debug_("selecting hyperslab with offsets: {}", utils::to_string(m_hyperslab_offsets));
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
        DEBUG_ASSERT_LT(iitem, m_nitem_global_max, "Item index exceeds global maximum");
        auto status = H5Sselect_none(m_filespace_handle);
        DEBUG_ONLY(status);
        DEBUG_ASSERT_FALSE(status, "HDF5 null selection failed");
    }
}

hdf5::NdListBase::~NdListBase() {
    H5Sclose(m_filespace_handle);
    H5Dclose(m_dataset_handle);
    H5Sclose(m_memspace_handle);
}

void hdf5::FileBase::check_is_hdf5(const std::string &name) {
    REQUIRE_TRUE(H5Fis_hdf5(name.c_str()), "Specified file is not HDF5 format");
}

hdf5::FileWriter::FileWriter(std::string name) : FileBase(H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, AccessPList())) {
    REQUIRE_GE(m_handle, 0,
                    "HDF5 file could not be opened for writing. It may be locked by another program");
}

hdf5::FileReader::FileReader(std::string name) : FileBase(
        (check_is_hdf5(name), H5Fopen(name.c_str(), H5F_ACC_RDONLY, AccessPList()))) {
    REQUIRE_GE(m_handle, 0, "HDF5 file could not be opened for reading.");
}
