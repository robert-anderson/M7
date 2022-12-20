//
// Created by rja on 06/12/22.
//

#ifndef M7_HDF5_BUFFERER_H
#define M7_HDF5_BUFFERER_H

#include <utility>

#include "Type.h"

/**
 *
 */
namespace hdf5 {
    namespace dataset {
        /**
         * nomenclature:
         *  type = an HDF5 native type (cannot be complex)
         *  elem = single number as stored in the C++ program (real or complex)
         *  h5_elem = single number as stored in the HDF5 archive (cannot be complex)
         *  item = a shaped array of elems
         *  h5_item = a shaped array of h5_elems
         */

        struct Format {
            /**
             * native type of the HDF5 dataset involved in the I/O
             */
            const Type m_h5_type;
            /**
             * layout of an item in C++ terms
             */
            const uintv_t m_shape;
            /**
             * layout of an item in HDF5 terms in hsize_t type
             */
            const v_t<hsize_t> m_h5_shape;
            /**
             * size in bytes of each item
             */
            const uint_t m_size;
            /**
             * names of each dimension of the item in the HDF5 metadata
             */
            const strv_t m_dim_names;

            Format(Type h5_type, uintv_t shape, strv_t dim_names, bool add_complex_dim);
        };

        /**
         * Adding an extra dimension to implement a concept of lists of multidimensional HDF5 arrays
         */
        struct ListFormat {
            /**
             * layout of a single item in the multidimensional list
             */
            Format m_item;
            /**
             * total number of items (locally held)
             */
            const uint_t m_nitem;
            /**
             * layout of all locally stored items in HDF5 terms in hsize_t type
             */
            const v_t<hsize_t> m_h5_shape;
            /**
             * names of each dimension in the HDF5 metadata
             */
            const strv_t m_dim_names;

            ListFormat(Format item_format, uint_t nitem);
        };

        struct DistListFormat {
            const ListFormat m_local;
            /**
             * sum of items across all ranks
             */
            const uint_t m_nitem;
            /**
             * index of start of items on this rank
             */
            const uint_t m_nitem_displ;
            /**
             * overall shape of distributed dataset with first element containing m_nitem
             */
            const v_t<hsize_t> m_h5_shape;

            DistListFormat(Format item_format, uint_t nitem);
        };

        typedef std::function<const void*(uint_t i, const DistListFormat& format, uint_t max_nitem_per_op)> save_fn;
        typedef std::function<void*(uint_t i, const DistListFormat& format, uint_t max_nitem_per_op)> load_fn;
    }
}

#endif //M7_HDF5_BUFFERER_H