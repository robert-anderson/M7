//
// Created by rja on 04/02/23.
//

#include "DatasetTransaction.h"


hdf5::DatasetTransaction::DatasetTransaction(hdf5::dataset::DistListFormat format) :
        m_format(std::move(format)), m_counts(m_format.m_h5_shape), m_offsets(m_format.m_local.ndim(), 0){
    // set number of transacted items to 0 - subclasses should set the initial element to the appropriate value
    // before each transaction (read/write call)
    m_counts[0] = 0;
}

hdf5::DatasetTransaction::~DatasetTransaction() {
    // free all HDF5 handles if they have been opened (set non-zero) by the subclasses
    if (m_plist) H5Pclose(m_plist);
    else logging::warn("property list was not opened in transaction");
    if(m_filespace) H5Sclose(m_filespace);
    else logging::warn("filespace was not opened in transaction");
    if (m_file_hyperslab) H5Sclose(m_file_hyperslab);
    else logging::warn("file hyperslab was not opened in transaction");
    if (m_mem_hyperslab) H5Sclose(m_mem_hyperslab);
    else logging::warn("application memory hyperslab was not opened in transaction");
    if (m_dataset) H5Dclose(m_dataset);
    else logging::warn("dataset was not opened in transaction");
}

hdf5::DatasetSaver::DatasetSaver(
        const hdf5::NodeWriter& nw, const str_t& name, hdf5::dataset::PartDistListFormat format) :
        DatasetTransaction(format){
    REQUIRE_FALSE_ALL(name.empty(), "HDF5 dataset must be given a name")
    m_filespace = H5Screate_simple(format.m_h5_shape.size(), format.m_h5_shape.data(), nullptr);

    // specify format of the dataset
    m_dataset = H5Dcreate(nw.m_handle, name.c_str(), format.m_local.m_item.m_type.m_handle, m_filespace,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    REQUIRE_GT_ALL(m_dataset, 0, "dataset creation failed");
    {
        // set dimension labels
        uint_t i = 0ul;
        for (auto& dim_name: format.m_local.m_dim_names) H5DSset_label(m_dataset, i++, dim_name.c_str());
    }
    H5Sclose(m_filespace);
    m_filespace = H5Dget_space(m_dataset);

    // datasets are always transacted in collective I/O mode
    m_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(m_plist, H5FD_MPIO_COLLECTIVE);

    const hsize_t ndim = format.m_h5_shape.size();
    // hyperslabs for the offsets and counts associated with the file and memory layout respectively
    m_file_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);
    m_mem_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);

    if (!m_format.m_nitem) logging::warn("Saving empty dataset \"{}\"", name);
}

bool hdf5::DatasetSaver::write(const void* src, uint_t nitem) {
    /*
     * only call H5Dwrite if there is any data across all MPI ranks
     */
    if (!m_format.m_nitem) return true;
    // the number of items transacted may not be larger than the number of items remaining
    REQUIRE_LE(m_nitem_done + nitem, m_format.m_local.m_nitem, "too many items");
    m_counts[0] = nitem;
    m_offsets[0] = 0;
    H5Sselect_hyperslab(m_mem_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
    m_offsets[0] = m_format.m_nitem_displ + m_nitem_done;
    REQUIRE_EQ(bool(src), bool(m_counts[0]), "count zero with non-null data or count non-zero with null data");
    H5Sselect_hyperslab(m_file_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
    auto status = H5Dwrite(m_dataset, m_format.m_local.m_item.m_type.m_handle, m_mem_hyperslab, m_file_hyperslab, m_plist, src);
    REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
    bool all_done = !src;
    all_done = mpi::all_land(all_done);
    // deplete number of remaining items by the number just transacted
    m_nitem_done += m_counts[0];
    return all_done;
}

void hdf5::DatasetSaver::save_dist_list(
    const NodeWriter& nw,
    const str_t& name,
    const void* src,
    Type type,
    bool is_complex,
    uintv_t item_shape,
    strv_t item_dim_names,
    uint_t nitem,
    std::list<Attr> attrs,
    uint_t max_nitem_per_op,
    str_t leading_dim_name)
{
    const dataset::ItemFormat item_format(type, std::move(item_shape), std::move(item_dim_names), is_complex);
    const dataset::PartDistListFormat format(item_format, nitem, leading_dim_name);
    DatasetSaver ds(nw, name, format);
    // reinterpret ptr as byte so arithmetic can be performed on it
    auto ptr = reinterpret_cast<const char*>(src);
    char all_done = false;
    while (!all_done) {
        const auto nitem_this_op = ds.nitem_next(max_nitem_per_op);
        all_done = ds.write(nitem_this_op ? ptr : nullptr, nitem_this_op);
        ptr += nitem_this_op * format.m_local.m_item.m_size;
    }
    ds.save_attrs(attrs);
}


void hdf5::DatasetSaver::save_attrs(const std::list<Attr>& attrs) {
    for (const auto& attr: attrs) attr.save(m_dataset);
}

bool hdf5::DatasetLoader::exists(hid_t parent, const str_t &name) {
    return H5Oexists_by_name(parent, name.c_str(), H5P_DEFAULT);
}

void hdf5::DatasetLoader::require_exists(hid_t parent, const str_t &name) {
    REQUIRE_TRUE_ALL(exists(parent, name), logging::format("Dataset \"{}\" does not exist", name));
}

uint_t hdf5::DatasetLoader::read_dataset_ndim(hid_t parent, const str_t &name) {
    require_exists(parent, name);
    auto dataset = H5Dopen1(parent, name.c_str());
    auto dataspace = H5Dget_space(dataset);
    auto rank = H5Sget_simple_extent_ndims(dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return rank;
}

uintv_t hdf5::DatasetLoader::read_shape(hid_t parent, const str_t &name) {
    require_exists(parent, name);
    auto ndim = read_dataset_ndim(parent, name);
    auto dataset = H5Dopen1(parent, name.c_str());
    REQUIRE_GT(dataset, 0, logging::format("no such dataset \"{}\"", name));
    auto dataspace = H5Dget_space(dataset);
    v_t<hsize_t> shape(ndim, 0ul);
    H5Sget_simple_extent_dims(dataspace, shape.data(), nullptr);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return convert::vector<uint_t>(shape);
}

hdf5::dataset::FullDistListFormat hdf5::DatasetLoader::read_full_format(hid_t parent, const str_t &name, bool this_rank) {
    REQUIRE_TRUE(H5Oexists_by_name(parent, name.c_str(), H5P_DEFAULT), "dataset is missing from HDF5 node");
    auto dataset = H5Dopen1(parent, name.c_str());
    const auto shape = read_shape(parent, name);
    const uint_t nitem = shape.front();
    const uintv_t item_shape(shape.cbegin() + 1, shape.cend());
    Type type(H5Dget_type(dataset));
    H5Dclose(dataset);
    return {{type, item_shape, {}, false}, this_rank ? nitem : 0ul};
}

hdf5::dataset::PartDistListFormat hdf5::DatasetLoader::read_part_format(hid_t parent, const str_t &name, bool this_rank) {
    const auto full_fmt = read_full_format(parent, name, this_rank);
    const auto &item = full_fmt.m_local.m_item;
    // reading is not required on all ranks, the boolean argument determines which ranks are involved in the read op
    auto reading_iranks = mpi::filter(this_rank);
    REQUIRE_FALSE(reading_iranks.empty(),
                  "cannot read on zero MPI ranks: this_rank must be true on at least one rank");
    uint_t nitem = 0ul;
    if (this_rank) {
        auto it = std::find(reading_iranks.cbegin(), reading_iranks.cend(), mpi::irank());
        DEBUG_ASSERT_EQ(it != reading_iranks.cend(), this_rank, "meant to be reading on this rank");
        const auto ibin = std::distance(reading_iranks.cbegin(), it);
        // share the reading workload evenly over the involved ranks
        nitem = integer::evenly_shared_count(full_fmt.m_nitem, ibin, reading_iranks.size());
    }
    return {{item.m_type, item.m_shape, item.m_dim_names, false}, nitem};
}

bool hdf5::DatasetLoader::valid_part_flag(bool part) {
    REQUIRE_TRUE(mpi::all_land(part) || mpi::all_land(!part),
                 "no mixing of partial and full reading is allowed, each rank must pass the same part value");
    return part;
}

hdf5::dataset::DistListFormat
hdf5::DatasetLoader::read_format(hid_t parent, const str_t &name, bool part, bool this_rank) {
    if (valid_part_flag(part)) return read_part_format(parent, name, this_rank);
    else return read_full_format(parent, name, this_rank);
}

hdf5::DatasetLoader::DatasetLoader(const hdf5::NodeReader &nr, const str_t &name, bool part, bool this_rank) :
        DatasetTransaction(read_format(nr.m_handle, name, part, this_rank)), m_name(name) {
    m_dataset = H5Dopen1(nr.m_handle, name.c_str());
    m_filespace = H5Dget_space(m_dataset);
    const auto ndim = m_format.m_h5_shape.size();

    // datasets are always transacted in collective I/O mode
    m_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(m_plist, H5FD_MPIO_COLLECTIVE);

    // hyperslabs for the offsets and counts associated with the file and memory layout respectively
    m_file_hyperslab = H5Screate_simple(ndim, m_format.m_h5_shape.data(), nullptr);
    m_mem_hyperslab = H5Screate_simple(ndim, m_format.m_h5_shape.data(), nullptr);

    if (!m_format.m_nitem) logging::warn("Loading empty dataset \"{}\"", name);
}

void hdf5::DatasetLoader::load_attrs(std::list<Attr> &attrs) const {
    attrs.clear();
    // load any attributes associated with the open dataset
    auto op = [](hid_t parent, const char *name, const H5A_info_t */*ainfo*/, void *attrs_cast) -> herr_t {
        auto attrs = reinterpret_cast<std::list<Attr> *>(attrs_cast);
        attrs->emplace_back(parent, name);
        return 0;
    };
    hsize_t idx = 0ul;
    // use creation order so that elementwise comparison between saved and loaded attr lists may be done
    H5Aiterate2(m_dataset, H5_INDEX_CRT_ORDER, H5_ITER_INC, &idx, op, reinterpret_cast<void *>(&attrs));
    REQUIRE_EQ(idx, attrs.size(), "number of attrs in list should match the final index");
}

bool hdf5::DatasetLoader::read(void *dst, uint_t nitem) {
    /*
     * only call H5Dread if there is any data across all MPI ranks
     */
    if (!m_format.m_nitem) return true;
    // every MPI rank must call read, since the I/O is done in collective mode
    mpi::barrier();
    // the number of items transacted may not be larger than the number of items remaining
    REQUIRE_LE(m_nitem_done + nitem, m_format.m_local.m_nitem, "too many items");
    m_counts[0] = nitem;
    m_offsets[0] = 0;
    H5Sselect_hyperslab(m_mem_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
    m_offsets[0] = m_format.m_nitem_displ + m_nitem_done;

    REQUIRE_EQ(bool(dst), bool(m_counts[0]),
               "nitem zero with non-null data or nitem non-zero with null data");
    H5Sselect_hyperslab(m_file_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
    auto status = H5Dread(m_dataset, m_format.m_local.m_item.m_type.m_handle, m_mem_hyperslab, m_file_hyperslab, m_plist, dst);
    REQUIRE_FALSE(status, "HDF5 Error on multidimensional load");

    bool all_done = !dst;
    all_done = mpi::all_land(all_done);
    // deplete number of remaining items by the number just transacted
    m_nitem_done += m_counts[0];
    return all_done;
}

hdf5::DatasetLoader::~DatasetLoader() {
    if (m_nitem_done && (m_nitem_done < m_format.m_local.m_nitem))
        logging::warn("Dataset \"{}\" was only partially loaded", m_name);
}

void hdf5::DatasetLoader::load_dist_list(
    const NodeReader& nr,
    const str_t& name,
    void* dst,
    uint_t size,
    hdf5::Type type,
    bool is_complex,
    bool part,
    bool this_rank,
    std::list<Attr>& attrs,
    uint_t max_nitem_per_op)
{
    DatasetLoader dl(nr, name, part, this_rank);
    REQUIRE_GE_ALL(size, dl.m_format.m_local.m_size,
                   logging::format("buffer is not large enough to read dataset \"{}\"", name));
    if (dl.m_format.m_local.m_size < size)
        logging::warn_("allocated buffer is over-sized for dataset \"{}\"", name);
    REQUIRE_TRUE_ALL(dl.m_format.m_local.m_item.m_type == type,
                   logging::format("native type of buffer does not agree with that of dataset \"{}\"", name));
    if (is_complex && dl.m_format.m_local.m_item.m_shape.back()!=2)
        logging::warn("buffer is complex but dataset \"{}\" does not have minor index of 2 (real/imag)", name);
    // reinterpret ptr as byte so arithmetic can be performed on it
    auto ptr = reinterpret_cast<char*>(dst);
    char all_done = false;
    while (!all_done) {
        const auto nitem_this_op = dl.nitem_next(max_nitem_per_op);
        all_done = dl.read(nitem_this_op ? ptr : nullptr, nitem_this_op);
        ptr += nitem_this_op * dl.m_format.m_local.m_item.m_size;
    }
    dl.load_attrs(attrs);
}