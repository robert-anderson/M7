//
// Created by rja on 22/12/22.
//

#include "M7_lib/util/Vector.h"
#include "NodeWriter.h"

void hdf5::NodeWriter::save_attr(const hdf5::Attr& attr) const {
    attr.save(m_handle);
}

void hdf5::NodeWriter::save_dataset(const str_t& name, dataset::save_fn fn, const dataset::PartDistListFormat& format,
                                    uint_t max_nitem_per_op, std::list<Attr> attrs) const {
    auto filespace = H5Screate_simple(format.m_h5_shape.size(), format.m_h5_shape.data(), nullptr);

    // specify format of the dataset
    auto dataset = H5Dcreate(m_handle, name.c_str(), format.m_local.m_item.m_type, filespace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    {
        // set dimension labels
        uint_t i = 0ul;
        for (auto& dim_name: format.m_local.m_dim_names) H5DSset_label(dataset, i++, dim_name.c_str());
    }
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
        const auto src = fn(format.m_local, max_nitem_per_op);
        REQUIRE_EQ(bool(src), bool(counts[0]), "count zero with non-null data or count non-zero with null data");
        H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts.data(), nullptr);
        auto status = H5Dwrite(dataset, format.m_local.m_item.m_type, mem_hyperslab, file_hyperslab, plist, src);
        REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
        all_done = !src;
        // only allow the loop to terminate if all ranks have yielded a null pointer
        all_done = mpi::all_land(all_done);
        // deplete number of remaining items by the number just transacted
        nitem_remaining -= counts[0];
    }

    // save any attributes before closing
    for (const auto& attr: attrs) attr.save(dataset);

    H5Pclose(plist);
    H5Sclose(filespace);
    H5Sclose(file_hyperslab);
    H5Sclose(mem_hyperslab);
    H5Dclose(dataset);
}


void hdf5::NodeWriter::save_dataset(const str_t& name, dataset::save_fn fn,
                                    const dataset::PartDistListFormat& format, std::list<Attr> attrs) const {
    save_dataset(name, fn, format, format.m_nitem, std::move(attrs));
}

void hdf5::NodeWriter::save_dataset(const str_t& name, dataset::save_fn fn,
                                    const dataset::PartDistListFormat& format) const {
    save_dataset(name, fn, format, format.m_nitem, {});
}