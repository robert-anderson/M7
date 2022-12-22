//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_NODE_H
#define M7_HDF5_NODE_H

#include "Attr.h"
#include "Dataset.h"
#include "IoManager.h"
#include "M7_lib/util/Pointer.h"

namespace hdf5 {

    struct Node {
        const hid_t m_handle;
        Node(hid_t handle);
        operator hid_t() const;
        bool attr_exists(const str_t& name) const;
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
         * @param lp
         *  determines loading behavior
         */
        template<typename T>
        void load_dataset(const str_t& name, T* dst, uint_t max_nitem_per_op, bool part, bool this_rank) const {
            auto fn = [&](uint_t i, const dataset::ListFormat& format, uint_t max_nitem_per_op) {
                auto begin = reinterpret_cast<char*>(dst);
                auto end = begin + format.m_size;
                auto ptr = begin + i * max_nitem_per_op * format.m_item.m_size;
                return ::ptr::in_range(ptr, begin, end) ? ptr : nullptr;
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

    struct NodeWriter : Node {
        NodeWriter(hid_t handle): Node(handle){}

        void save_attr(const str_t& name, const Attr& attr) const;

        template<typename T>
        void save_attr(const str_t& name, const T& v) const {
            save_attr(name, {v});
        }

        void save_dataset(const str_t& name, dataset::save_fn fn, const dataset::PartDistListFormat& format,
                          uint_t max_nitem_per_op) const;

        void save_dataset(const str_t& name, dataset::save_fn fn, const dataset::PartDistListFormat& format) const;

        /**
         * save data from a simple contiguous buffer pointed to by src into a HDF5 dataset.
         * @tparam T
         *  type of the source buffer
         * @param name
         *  name of the dataset
         * @param src
         *  source buffer
         * @param item_shape
         *  multidimensional shape of each item in the list
         * @param dim_names
         *  names of each dimension in the item
         * @param nitem
         *  number of items held on this rank
         * @param max_nitem_per_op
         *  maximum number of items in each transaction
         */
        template<typename T>
        void save_dataset(const str_t& name, const T* src, uintv_t item_shape,
                          strv_t dim_names, uint_t nitem, uint_t max_nitem_per_op) const {
            auto fn = [&](uint_t i, const dataset::ListFormat& format, uint_t max_nitem_per_op) {
                auto begin = reinterpret_cast<const char*>(src);
                auto end = begin + format.m_size;
                auto ptr = begin + i * max_nitem_per_op * format.m_item.m_size;
                return ::ptr::in_range(ptr, begin, end) ? ptr : nullptr;
            };
            return save_dataset(name, fn,
                    {{Type::make<T>(), item_shape, dim_names, dtype::is_complex<T>()}, nitem}, max_nitem_per_op);
        }

        template<typename T>
        void save_dataset(const str_t& name, const v_t<T>& src, uintv_t item_shape,
                          strv_t dim_names, uint_t max_nitem_per_op, bool this_rank) const {
            const auto nitem = this_rank ? src.size() / nd::nelement(item_shape) : 0ul;
            REQUIRE_FALSE(src.size() % nd::nelement(item_shape), "item shape inconsistent with size of data buffer");
            save_dataset(name, src.data(), item_shape, dim_names, nitem, max_nitem_per_op);
        }

        template<typename T>
        void save_dataset(const str_t& name, const v_t<T>& src, uintv_t item_shape, strv_t dim_names, bool this_rank) const {
            const auto nitem = this_rank ? src.size() / nd::nelement(item_shape) : 0ul;
            save_dataset(name, src.data(), item_shape, dim_names, nitem, nitem);
        }

        template<typename T>
        void save_dataset(const str_t& name, const v_t<T>& src, uint_t max_nitem_per_op, bool this_rank) const {
            save_dataset(name, src, {}, {}, max_nitem_per_op, this_rank);
        }

        template<typename T>
        void save_dataset(const str_t& name, const v_t<T>& src, bool this_rank) const {
            save_dataset(name, src, src.size(), this_rank);
        }

        template<typename T>
        void save_dataset(const str_t& name, T v, bool this_rank) const {
            v_t<T> src;
            src.push_back(v);
            save_dataset(name, src, this_rank);
        }
    };
}

#endif //M7_HDF5_NODE_H