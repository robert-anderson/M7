//
// Created by rja on 03/02/23.
//

#ifndef M7_DATASETLOADER_H
#define M7_DATASETLOADER_H

#include "DatasetFormat.h"
#include "NodeReader.h"
#include "M7_lib/util/Vector.h"

#if 0
namespace hdf5 {
    class DatasetLoader : public DatasetTransaction {
        static bool exists(hid_t parent, const str_t& name) {
            return H5Oexists_by_name(parent, name.c_str(), H5P_DEFAULT);
        }

        static void require_exists(hid_t parent, const str_t& name) {
            REQUIRE_TRUE_ALL(exists(parent, name),  logging::format("Dataset \"{}\" does not exist", name));
        }

        static uint_t read_dataset_ndim(hid_t parent, const str_t& name) {
            require_exists(parent, name);
            auto dataset = H5Dopen1(parent, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            auto rank = H5Sget_simple_extent_ndims(dataspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return rank;
        }

        static uintv_t read_shape(hid_t parent, const str_t& name) {
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

        static hdf5::dataset::FullDistListFormat read_full_format(hid_t parent, const str_t& name, bool this_rank) {
            REQUIRE_TRUE(H5Oexists_by_name(parent, name.c_str(), H5P_DEFAULT), "dataset is missing from HDF5 node");
            auto dataset = H5Dopen1(parent, name.c_str());
            const auto shape = read_shape(parent, name);
            const uint_t nitem = shape.front();
            const uintv_t item_shape(shape.cbegin()+1, shape.cend());
            Type type(H5Dget_type(dataset));
            H5Dclose(dataset);
            return {{type, item_shape, {}, false}, this_rank ? nitem : 0ul};
        }

        static hdf5::dataset::PartDistListFormat read_part_format(hid_t parent, const str_t& name, bool this_rank) {
            const auto full_fmt = read_full_format(parent, name, this_rank);
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
                nitem = integer::evenly_shared_count(full_fmt.m_nitem, ibin, reading_iranks.size());
            }
            return {{item.m_type, item.m_shape, item.m_dim_names, false}, nitem};
        }

        static bool valid_part_flag(bool part) {
            REQUIRE_TRUE(mpi::all_land(part) || mpi::all_land(!part),
                         "no mixing of partial and full reading is allowed, each rank must pass the same part value");
            return part;
        }

        /**
         * @param parent
         *  Node (group or file) open for reading
         * @param name
         *  name of the dataset whose sizes are queried
         * @param part
         *  true if this is a partial read, where each rank with this_rank=true has a portion of the data to read
         *  false if every rank with this_rank=true reads the entire dataset
         * @param this_rank
         *  true if this MPI rank participates in the load operation
         * @return
         *  pair of integers specifying the size of the locally-read data: nitem and item_size
         */
        static dataset::DistListFormat read_format(hid_t parent, const str_t& name, bool part, bool this_rank) {
            if (valid_part_flag(part)) return read_part_format(parent, name, this_rank);
            else return read_full_format(parent, name, this_rank);
        }


    public:

        DatasetLoader(const NodeReader& nr, const str_t& name, bool part, bool this_rank) :
                DatasetTransaction(read_format(nr.m_handle, name, part, this_rank)) {
            m_dataset = H5Dopen1(nr.m_handle, name.c_str());
            m_filespace = H5Dget_space(m_dataset);
            const auto ndim = m_format.m_h5_shape.size();

            // datasets are always transacted in collective I/O mode
            hid_t plist = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

            // initialize an array to specify the counts of data transacted. this only changes in the first element
            m_counts = vector::prepended(m_format.m_local.m_item.m_h5_shape, 0ul);
            // hyperslabs for the offsets and counts associated with the file and memory layout respectively
            m_file_hyperslab = H5Screate_simple(ndim, m_format.m_h5_shape.data(), nullptr);
            m_mem_hyperslab = H5Screate_simple(ndim, m_format.m_h5_shape.data(), nullptr);

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
                // get a pointer to which this portion of the dataset can be contiguously loaded
                const auto dst = prep_fn(format.m_local, max_nitem_per_op);
                REQUIRE_EQ(bool(dst), bool(counts[0]),
                           "nitem zero with non-null data or nitem non-zero with null data");
                H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, offsets.data(), nullptr, counts.data(), nullptr);
                auto status = H5Dread(dataset, format.m_local.m_item.m_type, mem_hyperslab, file_hyperslab, plist, dst);
                REQUIRE_FALSE(status, "HDF5 Error on multidimensional load");
                /*
                 * now that the buffer has been filled, let the target object be populated (this may be a null operation in the
                 * case that the contiguous buffer is directly part of the target object e.g. loading a std::vector)
                 */
                if (dst) fill_fn(dst, counts[0]);
                all_done = !dst;
                // only allow the loop to terminate if all ranks have yielded a null pointer
                all_done = mpi::all_land(all_done);
                // deplete number of remaining items by the number just transacted
                nitem_remaining -= counts[0];
            }

            // load any attributes associated with the open dataset
            {
                auto op = [](hid_t parent, const char* name, const H5A_info_t*/*ainfo*/, void* attrs_cast) -> herr_t {
                    auto attrs = reinterpret_cast<std::list<Attr>*>(attrs_cast);
                    attrs->emplace_back(parent, name);
                    return 0;
                };
                hsize_t idx = 0ul;
                // use creation order so that elementwise comparison between saved and loaded attr lists may be done
                H5Aiterate2(dataset, H5_INDEX_CRT_ORDER, H5_ITER_INC, &idx, op, reinterpret_cast<void*>(&attrs));
                REQUIRE_EQ(idx, attrs.size(), "number of attrs in list should match the final index");
            }

            H5Pclose(plist);
            H5Sclose(filespace);
            H5Sclose(file_hyperslab);
            H5Sclose(mem_hyperslab);
            H5Dclose(dataset);
        }
    };
}

#endif //M7_DATASETLOADER_H
#endif //M7_DATASETLOADER_H