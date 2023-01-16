//
// Created by rja on 22/12/22.
//

#ifndef M7_NODEREADER_H
#define M7_NODEREADER_H

#include "M7_lib/util/Pointer.h"
#include "Node.h"

namespace hdf5 {

    /**
     * Whereas saving an application object to one or more HDF5 datasets is a simple case of marshalling the data into a
     * contiguous buffer, loading datasets into application structures can be more complicated. This class manages this
     * process
     */
    struct DatasetLoader {
        const str_t m_name;
        const uint_t m_max_nitem_per_op;
        std::list<Attr> m_attrs;

        DatasetLoader(hid_t parent_handle, str_t name, const dataset::DistListFormat& format, uint_t max_nitem_per_op):
            m_name(std::move(name)) {


            REQUIRE_TRUE(H5Oexists_by_name(parent_handle, name.c_str(), )child_exists(name), "dataset is missing from HDF5 node");
            auto dataset = H5Dopen1(parent_handle, name.c_str());
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
                // get a pointer to which this portion of the dataset can be contiguously loaded
                const auto dst = prep_fn(format.m_local, max_nitem_per_op);
                REQUIRE_EQ(bool(dst), bool(counts[0]), "nitem zero with non-null data or nitem non-zero with null data");
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
                auto op = [](hid_t parent, const char* name, const H5A_info_t */*ainfo*/, void* attrs_cast) -> herr_t {
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

        void hdf5::NodeReader::load_dataset(const str_t& name, hdf5::dataset::load_prep_fn prep_fn,
                                            hdf5::dataset::load_fill_fn fill_fn,
                                            uint_t max_nitem_per_op, ) const {

        }

    };

    struct NodeReader : Node {

        NodeReader(hid_t handle) : Node(handle) {}

        bool child_exists(const str_t& name) const;

        uint_t first_existing_child(const strv_t& names) const;

        uint_t nchild() const;

        str_t child_name(uint_t ichild) const;

        int child_type(uint_t i) const;

        strv_t child_names(int type=-1) const;

        uint_t get_dataset_ndim(const str_t& name) const;

        uintv_t get_dataset_shape(const str_t& name) const;

        template<typename T>
        T load_attr(const str_t& name) const {
            T v;
            auto success = Attr(m_handle, name).parse(v);
            REQUIRE_TRUE(success, "HDF5 attribute load failed without default value");
            return v;
        }

        template<typename T>
        T load_attr(const str_t& name, T default_) const {
            T v;
            Attr(m_handle, name).parse(v, default_);
            return v;
        }

        dataset::FullDistListFormat get_full_dataset_format(const str_t& name, bool this_rank) const;

        dataset::PartDistListFormat get_part_dataset_format(const str_t& name, bool this_rank) const;

    private:
        static bool valid_part_flag(bool part);

        void get_local_sizes(const str_t& name, bool part, bool this_rank, uint_t& nitem, uint_t& item_size) const;

        template<typename T>
        uint_t local_nelement(const str_t& name, bool part, bool this_rank) const {
            uint_t nitem;
            uint_t item_size;
            get_local_sizes(name, part, this_rank, nitem, item_size);
            REQUIRE_FALSE(item_size % sizeof(T), "item size is invalid for the given type");
            return nitem * (item_size / sizeof(T));
        }

        void load_dataset(const str_t& name, dataset::load_prep_fn prep_fn, dataset::load_fill_fn fill_fn,
                          const dataset::DistListFormat& format, uint_t max_nitem_per_op, std::list<Attr>& attrs) const;

    public:
        void load_dataset(const str_t& name, dataset::load_prep_fn prep_fn, dataset::load_fill_fn fill_fn,
                          uint_t max_nitem_per_op, std::list<Attr>& attrs, bool part, bool this_rank) const;

        /**
         * load a dataset into a simple contiguous buffer pointed to by dst. it is obviously the responsibility of the
         * calling scope to ensure that dst points to memory of sufficient size to accommodate the data. The correct
         * size of the dst buffer can be determined by a preceding call to get_dataset_format
         * @tparam T
         *  type of the destination buffer
         * @param name
         *  name of the dataset
         * @param dst
         *  pointer to beginning of the destination buffer
         * @param max_nitem_per_op
         *  maximum number of items to transfer in a single operation
         * @param part
         *  true if the each included rank only loads part of the HDF5 dataset
         * @param this_rank
         *  true if this rank is active in the load operation
         */
        template<typename T>
        void load_dataset(const str_t& name, T* dst, uint_t max_nitem_per_op,
                          std::list<Attr>& attrs, bool part, bool this_rank) const {
            using namespace ::ptr;
            const auto begin = reinterpret_cast<buf_t*>(dst);
            auto ptr = begin;
            auto prep_fn = [&](const dataset::ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
                const auto tmp = in_range(ptr, begin, format.m_size);
                if (!tmp) return tmp;
                ptr = in_range(ptr + max_nitem_per_op * format.m_item.m_size, begin, format.m_size);
                return tmp;
            };
            // don't need a fill function here: prep_fn puts all data where it is supposed to be
            auto fill_fn = [](const buf_t*, uint_t){};
            return load_dataset(name, prep_fn, fill_fn, max_nitem_per_op, attrs, part, this_rank);
        }

        template<typename T>
        void load_dataset(const str_t& name, v_t<T>& dst, uint_t max_nitem_per_op,
                          std::list<Attr>& attrs, bool part, bool this_rank) const {
            dst.resize(local_nelement<T>(name, part, this_rank));
            load_dataset(name, dst.data(), max_nitem_per_op, attrs, part, this_rank);
        }

        template<typename T>
        void load_dataset(const str_t& name, v_t<T>& dst, std::list<Attr>& attrs, bool part, bool this_rank) const {
            dst.resize(local_nelement<T>(name, part, this_rank));
            load_dataset(name, dst.data(), dst.size(), attrs, part, this_rank);
        }

        template<typename T>
        T load_dataset(const str_t& name, std::list<Attr>& attrs, bool part, bool this_rank) const {
            T dst;
            load_dataset(name, dst, attrs, part, this_rank);
            return dst;
        }

        /*
         * below all the above methods are repeated with attributes discarded
         */

        void load_dataset(const str_t& name, dataset::load_prep_fn prep_fn, dataset::load_fill_fn fill_fn,
                          uint_t max_nitem_per_op, bool part, bool this_rank) const;

        template<typename T>
        void load_dataset(const str_t& name, T* dst, uint_t max_nitem_per_op, bool part, bool this_rank) const {
            std::list<Attr> attrs;
            load_dataset(name, dst, max_nitem_per_op, attrs, part, this_rank);
        }

        template<typename T>
        void load_dataset(const str_t& name, v_t<T>& dst, uint_t max_nitem_per_op, bool part, bool this_rank) const {
            std::list<Attr> attrs;
            load_dataset(name, dst, max_nitem_per_op, attrs, part, this_rank);
        }

        template<typename T>
        void load_dataset(const str_t& name, v_t<T>& dst, bool part, bool this_rank) const {
            std::list<Attr> attrs;
            load_dataset(name, dst, attrs, part, this_rank);
        }

        template<typename T>
        T load_dataset(const str_t& name, bool part, bool this_rank) const {
            std::list<Attr> attrs;
            return load_dataset<T>(name, attrs, part, this_rank);
        }
    };
}


#endif //M7_NODEREADER_H
