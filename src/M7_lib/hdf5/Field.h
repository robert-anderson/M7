//
// Created by rja on 23/12/22.
//

#ifndef M7_HDF5_FIELD_H
#define M7_HDF5_FIELD_H

#include <utility>

#include "M7_lib/field/FieldBase.h"
#include "Group.h"
#include "DatasetSaver.h"

namespace hdf5 {
    namespace field {
        static constexpr uint_t c_default_max_nitem_per_op = 16000000;
        static void save(const FieldBase& field, const hdf5::NodeWriter& nw, const str_t& name,
                         Type type, bool is_complex, uintv_t item_shape, strv_t item_dim_names, bool this_rank,
                         uint_t max_nitem_per_op = c_default_max_nitem_per_op, std::list<Attr> attrs = {}) {
            const dataset::ItemFormat item_format(type, std::move(item_shape), std::move(item_dim_names), is_complex);
            const dataset::PartDistListFormat format(item_format, this_rank ? field.m_row->m_table->nrecord() : 0ul, "row");
            REQUIRE_EQ(item_format.m_size, field.m_size, "item size implied by given shape inconsistent with field size");
            DatasetSaver ds(nw, name, format);

            // need not allocate a buffer longer than needed
            max_nitem_per_op = std::min(format.m_local.m_nitem, max_nitem_per_op);
            v_t<buf_t> buf(max_nitem_per_op * item_format.m_size);

            // items are given by records, which are rows below the high water mark that have not been freed
            uint_t nitem_found = 0ul;
            uint_t irow = 0ul;
            bool all_done = false;
            while (!all_done) {
                buf.clear();
                const auto next_nitem_found = std::min(nitem_found + max_nitem_per_op, format.m_local.m_nitem);
                if (next_nitem_found != nitem_found) {
                    auto buf_ptr = buf.data();
                    while (nitem_found != next_nitem_found) {
                        if (!field.m_row->m_table->is_freed(irow)) {
                            field.to_buffer(buf_ptr, irow);
                            buf_ptr += field.m_size;
                            ++nitem_found;
                        }
                        ++irow;
                    }
                }
                const void* src = (nitem_found==format.m_local.m_nitem) ? nullptr : buf.data();
                all_done = ds.write(src, std::min(ds.nitem_remaining(), max_nitem_per_op));
            }
            ds.save_attrs(attrs);
        }

        template<typename T>
        static void save(const FieldBase& field, const hdf5::NodeWriter& nw, const str_t& name,
                         uintv_t item_shape, strv_t item_dim_names, bool this_rank,
                         uint_t max_nitem_per_op = c_default_max_nitem_per_op, std::list<Attr> attrs = {}) {
            save(field, nw, name, Type::make<T>(), dtype::is_complex<T>(), item_shape,
                    item_dim_names, this_rank, max_nitem_per_op, attrs);
        }


        template<typename T>
        static void load(FieldBase& field, const hdf5::NodeReader& nr, const str_t& name,
                         uint_t max_nitem_per_op, std::list<Attr>& /*attrs*/, bool part, bool this_rank) {
            const auto local_format = part ? nr.get_part_dataset_format(name, this_rank).m_local :
                                      nr.get_full_dataset_format(name, this_rank).m_local;
            // unlike the save case, there is no need to consider empty rows in the table
            const uint_t nitem = local_format.m_nitem;

            v_t<buf_t> buf;
            uint_t iitem = 0ul;
            auto prep_fn = [&iitem, &nitem, &buf](const hdf5::dataset::ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
                if (iitem==nitem) return nullptr;
                if (buf.empty()) buf.resize(max_nitem_per_op * format.m_item.m_size);
                iitem = std::min(iitem + max_nitem_per_op, nitem);
                return buf.data();
            };

            // make sure the table has an adequate number of "in use" rows
            auto& table = field.m_row->m_table;
            // if there are not enough "in use" rows, make room
            if (table->nrow_in_use() < nitem) table->push_back(nitem-table->nrow_in_use());
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
