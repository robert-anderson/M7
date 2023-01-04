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
                         strv_t dim_names, uint_t max_nitem_per_op, const std::list<Attr>& attrs, bool this_rank) {
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
                         strv_t dim_names, const std::list<Attr>& attrs, bool this_rank) {
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

        template<typename T>
        static void load(FieldBase& field, const hdf5::NodeReader& nr, const str_t& name,
                         uint_t max_nitem_per_op, std::list<Attr>& /*attrs*/, bool part, bool this_rank) {
            const auto local_format = part ? nr.get_part_dataset_format(name, this_rank).m_local :
                                      nr.get_full_dataset_format(name, this_rank).m_local;
            // unlike the save case, there is no need to consider empty rows in the table
            const uint_t nitem = local_format.m_nitem;
            // make sure the table has an adequate number of "in use" rows
            auto& table = field.m_row->m_table;
            if (table->nrow_in_use() < nitem) table->push_back(nitem);
            v_t<buf_t> buf;
            uint_t iitem = 0ul;
            auto prep_fn = [&iitem, &nitem, &buf](const hdf5::dataset::ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
                if (iitem==nitem) return nullptr;
                if (buf.empty()) buf.resize(max_nitem_per_op * format.m_item.m_size);
                iitem = std::min(iitem + max_nitem_per_op, nitem);
                return buf.data();
            };

            uint_t irow = 0ul;
            auto fill_fn = [&irow, &field](const buf_t* src, uint_t nitem) -> void {
                const auto next_irow = irow + nitem;
                auto ptr = src;
                for (; irow < next_irow; ++irow) {
                    field.from_buffer(ptr, irow);
                    ptr +=  field.m_size;
                }
            };
            nr.load_dataset(name, prep_fn, fill_fn, max_nitem_per_op, part, this_rank);
        }
        /**
         * assume all items are loaded in one op
         */
        template<typename T>
        static void load(FieldBase& field, const hdf5::NodeReader& nr, const str_t& name,
                         std::list<Attr>& attrs, bool part, bool this_rank) {
            const auto local_format = part ? nr.get_part_dataset_format(name, this_rank).m_local :
                                      nr.get_full_dataset_format(name, this_rank).m_local;
            load<T>(field, nr, name, local_format.m_nitem, attrs, part, this_rank);
        }

        /**
         * assume no attributes
         */

        template<typename T>
        static void load(FieldBase& field, const hdf5::NodeReader& nr, const str_t& name,
                         uint_t max_nitem_per_op, bool part, bool this_rank) {
            std::list<Attr> attrs;
            load<T>(field, nr, name, max_nitem_per_op, attrs, part, this_rank);
        }

        /**
         * assume no attributes and all items are loaded in one op
         */
        template<typename T>
        static void load(FieldBase& field, const hdf5::NodeReader& nr, const str_t& name, bool part, bool this_rank) {
            std::list<Attr> attrs;
            load<T>(field, nr, name, attrs, part, this_rank);
        }
    }
}



#endif //M7_HDF5_FIELD_H
