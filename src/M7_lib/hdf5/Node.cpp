//
// Created by anderson on 27/06/2022.
//

#include "Node.h"
#include "M7_lib/util/Vector.h"

hdf5::Node::Node(hid_t handle) : m_handle(handle){}

hdf5::Node::operator hid_t() const {
    return m_handle;
}

bool hdf5::Node::attr_exists(const str_t& name) const {
    return H5Aexists(m_handle, name.c_str());
}

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

hdf5::dataset::DistListFormat hdf5::NodeReader::get_dataset_format(const str_t& name, LoadPolicy lp) const {
    char lp_is_all_load_all = (lp==AllLoadAll);
    REQUIRE_TRUE(mpi::all_land(lp_is_all_load_all) || mpi::all_land(!lp_is_all_load_all), "invalid LoadPolicy");

    REQUIRE_TRUE(child_exists(name), "dataset is missing from HDF5 node");
    auto dataset = H5Dopen1(m_handle, name.c_str());

    const auto shape = get_dataset_shape(name);
    const uint_t nitem_sum = shape.front();
    const uintv_t item_shape(shape.cbegin()+1, shape.cend());
    Type h5_type(H5Dget_type(dataset));
    H5Dclose(dataset);

    if (lp_is_all_load_all) return {{h5_type, item_shape, {}, false}, nitem_sum};

    const auto this_rank = (lp==PartialLoad);
    // reading is not required on all ranks, the boolean argument determines which ranks are involved in the read op
    auto reading_iranks = mpi::filter(this_rank);
    REQUIRE_FALSE(reading_iranks.empty(), "cannot read on zero MPI ranks: this_rank must be true on at least one rank");
    uint_t nitem = 0ul;
    if (this_rank) {
        auto it = std::find(reading_iranks.cbegin(), reading_iranks.cend(), mpi::irank());
        DEBUG_ASSERT_EQ(it != reading_iranks.cend(), this_rank, "meant to be reading on this rank");
        const auto ibin = std::distance(reading_iranks.cbegin(), it);
        // share the reading workload evenly over the involved ranks
        nitem = integer::evenly_shared_count(nitem_sum, ibin, reading_iranks.size());
    }
    return {{h5_type, item_shape, {}, false}, nitem};
}


void hdf5::NodeReader::load_dataset(const str_t& name, dataset::load_fn fn, uint_t max_nitem_per_op, LoadPolicy lp) const {
    REQUIRE_TRUE(child_exists(name), "dataset is missing from HDF5 node");
    auto dataset = H5Dopen1(m_handle, name.c_str());
    auto filespace = H5Dget_space(dataset);

    const auto format = get_dataset_format(name, lp);
    const auto ndim = format.m_h5_shape.size();

    // datasets are always transacted in collective I/O mode
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    // initialize an array to specify the counts of data transacted. this only changes in the first element
    auto counts = vector::prepended(format.m_local.m_item.m_h5_shape, 0ul);
    // hyperslabs for the offsets and counts associated with the file and memory layout respectively
    auto file_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);
    auto mem_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);

    if (lp==AllLoadAll)
        DEBUG_ASSERT_EQ(format.m_nitem, format.m_local.m_nitem,
                        "all local item counts should be equal to total item count");

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
        // do not add offset if the policy is AllLoadAll
        if (lp!=AllLoadAll) offsets[0] += format.m_nitem_displ;
        const auto dst = fn(iblock, format, max_nitem_per_op);
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

void hdf5::NodeWriter::save_attr(const str_t& name, const hdf5::Attr& attr) const {
    const auto& h5_shape = attr.m_format.m_h5_shape;
    auto dataspace = H5Screate_simple(h5_shape.size(), h5_shape.data(), nullptr);

    auto attr_handle = H5Acreate(m_handle, name.c_str(), attr.m_format.m_h5_type, dataspace,
                                 H5P_DEFAULT, H5P_DEFAULT);

    auto status = H5Awrite(attr_handle, attr.m_format.m_h5_type, attr.m_buf.data());
    DEBUG_ONLY(status);
    DEBUG_ASSERT_FALSE(status, "HDF5 attribute write failed");
    H5Aclose(attr_handle);
    H5Sclose(dataspace);
}

void hdf5::NodeWriter::save_dataset(const str_t& name, dataset::save_fn fn, const dataset::DistListFormat& format, uint_t max_nitem_per_op) const {
    auto filespace = H5Screate_simple(format.m_h5_shape.size(), format.m_h5_shape.data(), nullptr);

    // specify format of the dataset
    auto dataset = H5Dcreate(m_handle, name.c_str(), format.m_local.m_item.m_h5_type, filespace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    filespace = H5Dget_space(dataset);

    // datasets are always transacted in collective I/O mode
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    const hsize_t ndim = format.m_h5_shape.size();
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
        offsets[0] = std::min(iblock * max_nitem_per_op, format.m_local.m_nitem) + format.m_nitem_displ;
        const auto src = fn(iblock, format, max_nitem_per_op);
        REQUIRE_EQ(bool(src), bool(counts[0]), "count zero with non-null data or count non-zero with null data");
        H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts.data(), nullptr);
        auto status = H5Dwrite(dataset, format.m_local.m_item.m_h5_type, mem_hyperslab, file_hyperslab, plist, src);
        REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
        all_done = !src;
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


void hdf5::NodeWriter::save_dataset(const str_t& name, dataset::save_fn fn, const dataset::DistListFormat& format) const {
    save_dataset(name, fn, format, format.m_nitem);
}