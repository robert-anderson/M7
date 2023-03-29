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


        void load(
            FieldBase& field,
            const hdf5::NodeReader& nr,
            const str_t& name,
            uint_t max_nitem_per_op,
            std::list<Attr>& attrs,
            bool part,
            bool this_rank);

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
