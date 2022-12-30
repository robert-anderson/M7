//
// Created by rja on 23/12/22.
//

#ifndef M7_HDF5_FIELD_H
#define M7_HDF5_FIELD_H

#include "M7_lib/field/FieldBase.h"
#include "Group.h"

namespace hdf5 {
    namespace field {
        template<typename T>
        static void save(const FieldBase& field, const hdf5::NodeWriter& nw, const str_t& name, uintv_t item_shape,
                         strv_t dim_names, uint_t max_nitem_per_op, std::list<Attr> attrs, bool this_rank) {
            // items are given by records, which are rows below the high water mark that have not been freed
            const uint_t nitem = this_rank ? field.m_row->m_table->nrecord() : 0ul;
            uint_t nitem_found = 0ul;
            uint_t irow = 0ul;
            v_t<buf_t> buf;
            auto fn = [&](const hdf5::dataset::ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
                if (buf.empty()) buf.resize(max_nitem_per_op * format.m_item.m_size);
                buf.clear();
                const auto next_nitem_found = std::min(nitem_found + max_nitem_per_op, nitem);
                if (next_nitem_found != nitem_found) {
                    auto dst = buf.data();
                    while (nitem_found != next_nitem_found) {
                        if (!field.m_row->m_table->is_freed(irow)) {
                            field.to_buffer(dst, irow);
                            dst += field.m_size;
                            ++nitem_found;
                        }
                        ++irow;
                    }
                    return buf.data();
                }
                return nullptr;
            };
            REQUIRE_EQ(nd::nelement(item_shape)*sizeof(T), field.m_size, "");
            const hdf5::dataset::ItemFormat item_format(hdf5::Type::make<T>(), item_shape, dim_names, false);
            nw.save_dataset(name, fn, {item_format, nitem}, max_nitem_per_op, std::move(attrs));
        }

        /**
         * assume all items are saved in one op
         */
        template<typename T>
        static void save(const FieldBase& field, const hdf5::NodeWriter& nw, const str_t& name, uintv_t item_shape,
                         strv_t dim_names, std::list<Attr> attrs, bool this_rank) {
            const uint_t nitem = this_rank ? field.m_row->m_table->nrecord() : 0ul;
            save<T>(field, nw, name, item_shape, dim_names, nitem, attrs, this_rank);
        }
        /**
         * assume no attributes
         */
        template<typename T>
        static void save(const FieldBase& field, const hdf5::NodeWriter& nw, const str_t& name, uintv_t item_shape,
                         strv_t dim_names, uint_t max_nitem_per_op, bool this_rank) {
            save<T>(field, nw, name, item_shape, dim_names, max_nitem_per_op, {}, this_rank);
        }

        /**
         * assume no attributes and all items are saved in one op
         */
        template<typename T>
        static void save(const FieldBase& field, const hdf5::NodeWriter& nw, const str_t& name, uintv_t item_shape,
                         strv_t dim_names, bool this_rank) {
            save<T>(field, nw, name, item_shape, dim_names, {}, this_rank);
        }

#if 0
        template<typename T, typename postload_fn_t>
        static void load(FieldBase& field, const hdf5::NodeReader& nr, postload_fn_t postload_fn,
                         const str_t& name, bool this_rank) {
            functor::assert_prototype<void(T* /*data*/, uint_t /*irow*/)>(postload_fn);
            // items are given by records, which are rows below the high water mark that have not been freed
            const uint_t nitem = this_rank ? field.m_row->m_table->nrecord() : 0ul;
            uint_t nitem_found = 0ul;
            uint_t irow = 0ul;
            v_t<buf_t> buf;
            auto fn = [&](const hdf5::dataset::ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
                if (buf.empty()) buf.resize(max_nitem_per_op * format.m_item.m_size);
                buf.clear();
                const auto next_nitem_found = std::min(nitem_found + max_nitem_per_op, nitem);
                if (next_nitem_found != nitem_found) {
                    auto dst = buf.data();
                    while (nitem_found != next_nitem_found) {
                        if (!field.m_row->m_table->is_freed(irow)) {
                            field.to_buffer(dst, irow);
                            dst += field.m_size;
                            ++nitem_found;
                        }
                        ++irow;
                    }
                    presave_fn(reinterpret_cast<T*>(buf.data()), irow);
                    return buf.data();
                }
                return nullptr;
            };
            hdf5::dataset::ItemFormat item_format(hdf5::Type::make<T>(), {field.m_size / sizeof(T)}, {"bytes"}, false);
            nr.load_dataset(name, fn, {item_format, nitem});
        }
#endif
    }
}



#endif //M7_HDF5_FIELD_H
