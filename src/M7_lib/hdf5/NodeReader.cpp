//
// Created by rja on 22/12/22.
//
#include "M7_lib/util/Vector.h"
#include "NodeReader.h"

hdf5::Attr hdf5::NodeReader::load_raw_attr(const str_t& name) const {
    if (!attr_exists(name)) return {{}, {}};
    auto attr_handle = H5Aopen(m_handle, name.c_str(), H5P_DEFAULT);
    auto dataspace = H5Aget_space(attr_handle);
    auto ndim = H5Sget_simple_extent_ndims(dataspace);
    v_t<hsize_t> shape(ndim);
    H5Sget_simple_extent_dims(dataspace, shape.data(), nullptr);
    Type h5_type(H5Aget_type(attr_handle));
    dataset::Format format(h5_type, convert::vector<uint_t>(shape), {}, false);
    v_t<buf_t> buf(format.m_size);
    auto status = H5Aread(attr_handle, h5_type, buf.data());
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute write failed");
    H5Aclose(attr_handle);
    H5Sclose(dataspace);
    return {buf, format};
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
    v_t<hsize_t> shape(ndim, 0ul);
    H5Sget_simple_extent_dims(dataspace, shape.data(), nullptr);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return convert::vector<uint_t>(shape);
}

hdf5::dataset::FullDistListFormat hdf5::NodeReader::get_full_dataset_format(const str_t& name, bool this_rank) const {
    REQUIRE_TRUE(child_exists(name), "dataset is missing from HDF5 node");
    auto dataset = H5Dopen1(m_handle, name.c_str());

    const auto shape = get_dataset_shape(name);
    const uint_t nitem = shape.front();
    const uintv_t item_shape(shape.cbegin()+1, shape.cend());
    Type h5_type(H5Dget_type(dataset));
    H5Dclose(dataset);
    return {{h5_type, item_shape, {}, false}, this_rank ? nitem : 0ul};
}

hdf5::dataset::PartDistListFormat hdf5::NodeReader::get_part_dataset_format(const str_t& name, bool this_rank) const {
    const auto full_fmt = get_full_dataset_format(name, this_rank);
    const auto& item = full_fmt.m_local.m_item;
    // reading is not required on all ranks, the boolean argument determines which ranks are involved in the read op
    auto reading_iranks = mpi::filter(this_rank);
    REQUIRE_FALSE(reading_iranks.empty(), "cannot read on zero MPI ranks: this_rank must be true on at least one rank");
    uint_t nitem = 0ul;
    if (this_rank) {
        auto it = std::find(reading_iranks.cbegin(), reading_iranks.cend(), mpi::irank());
        DEBUG_ASSERT_EQ(it != reading_iranks.cend(), this_rank, "meant to be reading on this rank");
        const auto ibin = std::distance(reading_iranks.cbegin(), it);
        // share the reading workload evenly over the involved ranks
        nitem = integer::evenly_shared_count(nitem, ibin, reading_iranks.size());
    }
    return {{item.m_h5_type, item.m_shape, item.m_dim_names, false}, nitem};
}

void hdf5::NodeReader::load_dataset(const str_t& name, hdf5::dataset::load_fn fn,
                                    const hdf5::dataset::DistListFormat& format, uint_t max_nitem_per_op) const {
    REQUIRE_TRUE(child_exists(name), "dataset is missing from HDF5 node");
    auto dataset = H5Dopen1(m_handle, name.c_str());
    auto filespace = H5Dget_space(dataset);
    const auto ndim = format.m_h5_shape.size();

    // datasets are always transacted in collective I/O mode
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    // initialize an array to specify the counts of data transacted. this only changes in the first element
    auto counts = vector::prepended(format.m_local.m_item.m_h5_shape, 0ul);
    // hyperslabs for the offsets and counts associated with the file and memory layout respectively
    auto file_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);
    auto mem_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);

    // number of items yet to be transacted
    uint_t nitem_remaining = format.m_local.m_nitem;
    // first element of offset is incremented at each iteration
    v_t<hsize_t> offsets(ndim);
    char all_done = false;
    for (uint_t iblock = 0ul; !all_done; ++iblock) {
        // the number of items transacted may not be larger than the cutoff
        counts[0] = std::min(max_nitem_per_op, nitem_remaining);
        offsets[0] = 0;
        H5Sselect_hyperslab(mem_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts.data(), nullptr);
        offsets[0] = std::min(iblock * max_nitem_per_op, format.m_local.m_nitem);
        offsets[0] += format.m_nitem_displ;
        const auto dst = fn(iblock, format.m_local, max_nitem_per_op);
        REQUIRE_EQ(bool(dst), bool(counts[0]), "nitem zero with non-null data or nitem non-zero with null data");
        H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts.data(), nullptr);
        auto status = H5Dread(dataset, format.m_local.m_item.m_h5_type, mem_hyperslab, file_hyperslab, plist, dst);
        REQUIRE_FALSE(status, "HDF5 Error on multidimensional load");
        all_done = !dst;
        // only allow the loop to terminate if all ranks have yielded a null pointer
        all_done = mpi::all_land(all_done);
        // deplete number of remaining items by the number just transacted
        nitem_remaining -= counts[0];
    }

    H5Pclose(plist);
    H5Sclose(filespace);
    H5Sclose(file_hyperslab);
    H5Sclose(mem_hyperslab);
    H5Dclose(dataset);
}