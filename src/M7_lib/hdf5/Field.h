//
// Created by rja on 23/12/22.
//

#ifndef M7_HDF5_FIELD_H
#define M7_HDF5_FIELD_H

#include <utility>

#include "M7_lib/field/FieldBase.h"
#include "Node.h"
#include "DatasetTransaction.h"

namespace hdf5 {
    namespace field {
        void save(
            const FieldBase& field,
            const hdf5::NodeWriter& nw,
            const str_t& name,
            Type type,
            bool is_complex,
            uintv_t item_shape,
            strv_t item_dim_names,
            bool this_rank,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op,
            std::list<Attr> attrs = {});

        template<typename T>
        static void save(
            const FieldBase& field,
            const hdf5::NodeWriter& nw,
            const str_t& name,
            uintv_t item_shape,
            strv_t item_dim_names,
            bool this_rank,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op,
            std::list<Attr> attrs = {})
        {
            save(field, nw, name, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names,
                 this_rank, max_nitem_per_op, attrs);
        }


        static void load(
            FieldBase& field,
            const hdf5::NodeReader& nr,
            const str_t& name,
            uint_t max_nitem_per_op,
            std::list<Attr>& attrs,
            bool part,
            bool this_rank)
        {
            DatasetLoader dl(nr, name, part, this_rank);
            if (field.m_row->m_table->nrow_in_use() < dl.nitem_remaining())
                field.m_row->m_table->push_back(dl.nitem_remaining());

            // need not allocate a buffer longer than needed
            max_nitem_per_op = std::min(dl.m_format.m_local.m_nitem, max_nitem_per_op);
            v_t<buf_t> buf(max_nitem_per_op * dl.m_format.m_local.m_item.m_size);

            // items are given by records, which are rows below the high water mark that have not been freed
            uint_t nitem_found = 0ul;
            uint_t irow = 0ul;
            bool all_done = false;
            while (!all_done) {
                buf.clear();
                const auto nitem_to_find = dl.nitem_next(max_nitem_per_op);
                const auto next_nitem_found = nitem_found + nitem_to_find;
                all_done = dl.read(nitem_to_find ? buf.data() : nullptr, nitem_to_find);
                if (next_nitem_found != nitem_found) {
                    auto buf_ptr = buf.data();
                    while (nitem_found != next_nitem_found) {
                        if (!field.m_row->m_table->is_freed(irow)) {
                            field.from_buffer(buf_ptr, irow);
                            buf_ptr += field.m_size;
                            ++nitem_found;
                        }
                        ++irow;
                    }
                }
            }
            dl.load_attrs(attrs);
        }

        /**
         * assume all items are loaded in one op
         */
        static void load(
            FieldBase& field,
            const hdf5::NodeReader& nr,
            const str_t& name,
            std::list<Attr>& attrs,
            bool part,
            bool this_rank)
        {
            const auto nitem = DatasetLoader::read_format(nr.m_id, name, part, this_rank).m_local.m_nitem;
            load(field, nr, name, nitem, attrs, part, this_rank);
        }

        /**
         * assume no attributes
         */
        static void load(
            FieldBase& field,
            const hdf5::NodeReader& nr,
            const str_t& name,
            uint_t max_nitem_per_op,
            bool part,
            bool this_rank)
        {
            std::list<Attr> attrs;
            load(field, nr, name, max_nitem_per_op, attrs, part, this_rank);
        }

        /**
         * assume no attributes and all items are loaded in one op
         */
        static void load(FieldBase& field, const hdf5::NodeReader& nr, const str_t& name, bool part, bool this_rank) {
            std::list<Attr> attrs;
            load(field, nr, name, attrs, part, this_rank);
        }
    }
}



#endif //M7_HDF5_FIELD_H
