//
// Created by anderson on 27/06/2022.
//

#include "Node.h"

hdf5::Node::Node(hid_t handle) : m_handle(handle){}

hdf5::Node::operator hid_t() const {
    return m_handle;
}

bool hdf5::Node::attr_exists(const str_t& name) const {
    return H5Aexists(m_handle, name.c_str());
}

bool hdf5::NodeReader::child_exists(const str_t &name) const {
    for (uint_t i=0ul; i<nchild(); ++i) if (name==child_name(i)) return true;
    return false;
}

uint_t hdf5::NodeReader::first_existing_child(const strv_t &names) const {
    for (uint_t i=0ul; i<names.size(); ++i) if (child_exists(names[i])) return i;
    return ~0ul;
}

uint_t hdf5::NodeReader::nchild() const {
    hsize_t n;
    auto status = H5Gget_num_objs(m_handle, &n);
    REQUIRE_TRUE(!status, "could not get number of objects within HDF5 group");
    return n;
}

str_t hdf5::NodeReader::child_name(uint_t ichild) const {
    uint_t size = H5Lget_name_by_idx(m_handle, ".", H5_INDEX_NAME,
                                     H5_ITER_INC, ichild, nullptr, 0, H5P_DEFAULT);
    str_t name(size, 0);
    auto name_ptr = const_cast<char*>(name.c_str());
    uint_t size_chk = H5Lget_name_by_idx(m_handle, ".", H5_INDEX_NAME,
                                         H5_ITER_INC, ichild, name_ptr, size+1, H5P_DEFAULT);
    REQUIRE_EQ(size, size_chk, "inconsistent number of chars in name");
    return name;
}

int hdf5::NodeReader::child_type(uint_t i) const {
    return H5Gget_objtype_by_idx(m_handle, i);
}

strv_t hdf5::NodeReader::child_names(int type) const {
    strv_t names;
    auto n = nchild();
    names.reserve(n);
    for (uint_t i=0ul; i<n; ++i) {
        if (type<0 || child_type(i)==type) names.push_back(child_name(i));
    }
    return names;
}

uint_t hdf5::NodeReader::get_dataset_ndim(const str_t& name) const {
    auto status = H5Gget_objinfo(m_handle, name.c_str(), 0, nullptr);
    REQUIRE_TRUE(!status, "Dataset \"" + name + "\" does not exist");
    auto dataset = H5Dopen1(m_handle, name.c_str());
    auto dataspace = H5Dget_space(dataset);
    auto rank = H5Sget_simple_extent_ndims(dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return rank;
}

uintv_t hdf5::NodeReader::get_dataset_shape(const str_t& name) const {
    auto ndim = get_dataset_ndim(name);
    auto dataset = H5Dopen1(m_handle, name.c_str());
    REQUIRE_GT_ALL(dataset, 0, logging::format("no such dataset \"{}\"", name));
    auto dataspace = H5Dget_space(dataset);
    v_t<hsize_t> dims(ndim, 0ul);
    H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    uintv_t out;
    out.reserve(dims.size());
    for (const auto &i: dims) out.push_back(i);
    return out;
}

void hdf5::NodeWriter::save_dataset(const str_t& name, const hdf5::io_manager::SaveManager& manager) const {
    const auto& format = manager.m_format;
    auto filespace = H5Screate_simple(format.m_h5_shape_sum.size(), format.m_h5_shape_sum.data(), nullptr);

    /*
     * Create the dataset with default properties and close filespace.
     */
    auto dataset_handle = H5Dcreate(m_handle, name.c_str(), format.m_h5_type, filespace,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//                m_dataset_handle = H5Dopen1(m_parent_handle, name.c_str());

    H5Sclose(filespace);

    filespace = H5Dget_space(dataset_handle);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    const hsize_t ndim = format.m_h5_shape.size();
    /*
     * we require three different hyperslabs:
     *  full: for writing a block of the maximum size
     *  part: for writing a block of the final non-zero size
     *  null: for writing nothing
     */
    const auto full_counts = vector::prepended(format.m_h5_item_shape, manager.m_max_nitem_per_write);
    const auto part_counts = vector::prepended(format.m_h5_item_shape, format.m_nitem % manager.m_max_nitem_per_write);
    const v_t<hsize_t> null_counts(ndim, 0);
    auto file_hyperslab = H5Screate_simple(ndim, manager.m_format.m_h5_shape_sum.data(), nullptr);
    auto mem_hyperslab = H5Screate_simple(ndim, manager.m_format.m_h5_shape_sum.data(), nullptr);

    uint_t nitem_remaining = format.m_nitem;
    v_t<hsize_t> offsets(ndim);
    char all_done = false;
    for (uint_t iblock=0ul; !all_done; ++iblock) {
        offsets[0] = 0;
        const hsize_t* counts = (nitem_remaining) ? (nitem_remaining < manager.m_max_nitem_per_write ? part_counts.data() : full_counts.data()) : null_counts.data();
        H5Sselect_hyperslab(mem_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts, nullptr);
        offsets[0] = std::min(iblock*manager.m_max_nitem_per_write, format.m_nitem) + format.m_nitem_offsets[mpi::irank()];
        const auto src = manager.get(iblock);
        REQUIRE_EQ(bool(src), bool(counts[0]), "count zero with non-null data or count non-zero with null data");
        H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts, nullptr);
        auto status = H5Dwrite(dataset_handle, format.m_h5_type, mem_hyperslab, file_hyperslab, plist, src);
        REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
        all_done = !src;
        all_done = mpi::all_land(all_done);
        nitem_remaining = (nitem_remaining > manager.m_max_nitem_per_write) ? nitem_remaining - manager.m_max_nitem_per_write : 0ul;
    }

    H5Sclose(filespace);
    H5Sclose(file_hyperslab);
    H5Sclose(mem_hyperslab);
    H5Dclose(dataset_handle);
}
