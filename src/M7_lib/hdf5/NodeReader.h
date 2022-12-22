//
// Created by rja on 22/12/22.
//

#ifndef M7_NODEREADER_H
#define M7_NODEREADER_H

#include "Node.h"

namespace hdf5 {

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

    private:
        Attr load_raw_attr(const str_t& name) const;

    public:

        template<typename T>
        T load_attr(const str_t& name) const {
            T v;
            auto success = load_raw_attr(name).parse(v);
            REQUIRE_TRUE(success, "HDF5 attribute load failed without default value");
            return v;
        }

        template<typename T>
        T load_attr(const str_t& name, T default_) const {
            T v;
            if (!load_raw_attr(name).parse(v)) v = default_;
            return v;
        }

        dataset::FullDistListFormat get_full_dataset_format(const str_t& name, bool this_rank) const;

        dataset::PartDistListFormat get_part_dataset_format(const str_t& name, bool this_rank) const;

    private:
        static bool valid_part_flag(bool part) {
            REQUIRE_TRUE(mpi::all_land(part) || mpi::all_land(!part),
                         "no mixing of partial and full reading is allowed, each rank must pass the same part value");
            return part;
        }

        void get_local_sizes(const str_t& name, bool part, bool this_rank, uint_t& nitem, uint_t& item_size) const {
            if (valid_part_flag(part)) {
                const auto format = get_part_dataset_format(name, this_rank);
                nitem = format.m_local.m_nitem;
                item_size = format.m_local.m_item.m_size;
            } else {
                const auto format = get_full_dataset_format(name, this_rank);
                nitem = format.m_local.m_nitem;
                item_size = format.m_local.m_item.m_size;
            }
        }

        template<typename T>
        uint_t local_nelement(const str_t& name, bool part, bool this_rank) const {
            uint_t nitem;
            uint_t item_size;
            get_local_sizes(name, part, this_rank, nitem, item_size);
            REQUIRE_FALSE(item_size % sizeof(T), "item size is invalid for the given type");
            return nitem * (item_size / sizeof(T));
        }

        void load_dataset(const str_t& name, dataset::load_fn fn, const dataset::DistListFormat& format, uint_t max_nitem_per_op) const;

    public:
        void load_dataset(const str_t& name, dataset::load_fn fn, uint_t max_nitem_per_op, bool part, bool this_rank) const {
            if (valid_part_flag(part)) {
                const auto format = get_part_dataset_format(name, this_rank);
                // dispatch the private method
                load_dataset(name, fn, format, max_nitem_per_op);
            } else {
                const auto format = get_full_dataset_format(name, this_rank);
                // dispatch the private method
                load_dataset(name, fn, format, max_nitem_per_op);
            }
        }

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
        void load_dataset(const str_t& name, T* dst, uint_t max_nitem_per_op, bool part, bool this_rank) const {
            using namespace ptr;
            const auto begin = reinterpret_cast<char*>(dst);
            auto ptr = begin;
            auto fn = [&](const dataset::ListFormat& format, uint_t max_nitem_per_op) {
                const auto tmp = in_range(ptr, begin, format.m_size);
                if (!tmp) return tmp;
                ptr = in_range(ptr + max_nitem_per_op * format.m_item.m_size, begin, format.m_size);
                return tmp;
            };
            return load_dataset(name, fn, max_nitem_per_op, part, this_rank);
        }

        template<typename T>
        void load_dataset(const str_t& name, v_t<T>& dst, uint_t max_nitem_per_op, bool part, bool this_rank) const {
            dst.resize(local_nelement<T>(name, part, this_rank));
            load_dataset(name, dst.data(), max_nitem_per_op, part, this_rank);
        }

        template<typename T>
        void load_dataset(const str_t& name, v_t<T>& dst, bool part, bool this_rank) const {
            dst.resize(local_nelement<T>(name, part, this_rank));
            load_dataset(name, dst.data(), dst.size(), part, this_rank);
        }

        template<typename T>
        T load_dataset(const str_t& name, bool part, bool this_rank) const {
            T dst;
            load_dataset(name, dst, part, this_rank);
            return dst;
        }
    };
}


#endif //M7_NODEREADER_H
