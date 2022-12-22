//
// Created by rja on 22/12/22.
//

#ifndef M7_NODEWRITER_H
#define M7_NODEWRITER_H

#include "Node.h"

namespace hdf5 {

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


#endif //M7_NODEWRITER_H
