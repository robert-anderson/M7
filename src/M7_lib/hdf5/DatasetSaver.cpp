//
// Created by rja on 30/01/23.
//

#include "DatasetSaver.h"

hdf5::DatasetSaver::DatasetSaver(
        const hdf5::NodeWriter& nw, const str_t& name, hdf5::dataset::PartDistListFormat format) :
        m_format(std::move(format)){
    REQUIRE_FALSE_ALL(name.empty(), "HDF5 dataset must be given a name")
    m_filespace = H5Screate_simple(format.m_h5_shape.size(), format.m_h5_shape.data(), nullptr);

    // specify format of the dataset
    m_dataset = H5Dcreate(nw.m_handle, name.c_str(), format.m_local.m_item.m_type, m_filespace,
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
    // initialize an array to specify the counts of data transacted. this only changes in the first element
    m_counts = vector::prepended(format.m_local.m_item.m_h5_shape, 0ul);
    // hyperslabs for the offsets and counts associated with the file and memory layout respectively
    m_file_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);
    m_mem_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);

    // first element of offset is incremented at each iteration
    m_offsets = v_t<hsize_t>(ndim);

    if (!m_format.m_nitem) logging::warn("Saving empty dataset \"{}\"", name);
}

void hdf5::DatasetSaver::save_dist_list(
        const hdf5::NodeWriter& nw, const str_t& name, const void* src, hdf5::Type type, bool is_complex,
        uintv_t item_shape, strv_t item_dim_names, uint_t nitem, std::list<Attr> attrs, str_t leading_dim_name) {
    const dataset::ItemFormat item_format(type, std::move(item_shape), std::move(item_dim_names), is_complex);
    const dataset::PartDistListFormat format(item_format, nitem, std::move(leading_dim_name));
    DatasetSaver ds(nw, name, format);
    ds.write(src, nitem);
    ds.save_attrs(attrs);
}

bool hdf5::DatasetSaver::write(const void* src, uint_t nitem) {
    /*
     * only call H5Dwrite if there is any data across all MPI ranks
     */
    if (!m_format.m_nitem) return true;
    // the number of items transacted may not be larger than the number of items remaining
    REQUIRE_LE(m_nitem_saved + nitem, m_format.m_local.m_nitem, "too many items");
    m_counts[0] = nitem;
    m_offsets[0] = 0;
    H5Sselect_hyperslab(m_mem_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
    m_offsets[0] = m_format.m_nitem_displ + m_nitem_saved;
    REQUIRE_EQ(bool(src), bool(m_counts[0]), "count zero with non-null data or count non-zero with null data");
    H5Sselect_hyperslab(m_file_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
    auto status = H5Dwrite(m_dataset, m_format.m_local.m_item.m_type, m_mem_hyperslab, m_file_hyperslab, m_plist, src);
    REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
    bool all_done = !src;
    all_done = mpi::all_land(all_done);
    // deplete number of remaining items by the number just transacted
    m_nitem_saved += m_counts[0];
    return all_done;
}

void hdf5::DatasetSaver::save_attrs(const std::list<Attr>& attrs) {
    for (const auto& attr: attrs) attr.save(m_dataset);
}

hdf5::DatasetSaver::~DatasetSaver() {
    // free all HDF5 handles
    H5Pclose(m_plist);
    H5Sclose(m_filespace);
    H5Sclose(m_file_hyperslab);
    H5Sclose(m_mem_hyperslab);
    H5Dclose(m_dataset);
}
