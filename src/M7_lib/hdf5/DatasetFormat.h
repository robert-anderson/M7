//
// Created by rja on 06/12/22.
//

#ifndef M7_DATASET_FORMAT_H
#define M7_DATASET_FORMAT_H

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

        struct ItemFormat {
            /**
             * native type of the HDF5 dataset involved in the I/O
             */
            const Type m_type;
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

            ItemFormat(Type type, uintv_t shape, strv_t dim_names, bool add_complex_dim);

            ItemFormat() : ItemFormat(Type(), {}, {}, false){}

            bool operator==(const ItemFormat& other) const;

            bool operator!=(const ItemFormat& other) const;
        };

        /**
         * Adding an extra dimension to implement a concept of lists of multidimensional HDF5 arrays
         */
        struct ListFormat {
            /**
             * layout of a single item in the multidimensional list
             */
            ItemFormat m_item;
            /**
             * total number of items (locally held)
             */
            const uint_t m_nitem;
            /**
             * total number of bytes (locally held)
             */
            const uint_t m_size;
            /**
             * layout of all locally stored items in HDF5 terms in hsize_t type
             */
            const v_t<hsize_t> m_h5_shape;
            /**
             * names of each dimension in the HDF5 metadata
             */
            const strv_t m_dim_names;

            /**
             * @return
             *  number of dimensions including the list dimension (dimensionality of item format + 1)
             */
            uint_t ndim() const;

            ListFormat(ItemFormat item_format, uint_t nitem, str_t leading_dim_name="item");
        };
        /**
         * the save-load functions only need access to the list format, not the details of the distribution
         *
         * save prepares a contiguous buffer representing the data to be saved. the caller of save_fn only needs to
         * perform the write operation
         *
         * load requires two functions:
         *  - one to prepare the pointer to the position in the buffer into which the next hyperslab from the HDF5
         *    dataset should be loaded
         *  - another to fill this buffered data into the target object
         */
        typedef std::function<const buf_t*(const ListFormat& format, uint_t max_nitem_per_op)> save_fn;
        typedef std::function<buf_t*(const ListFormat& format, uint_t max_nitem_per_op)> load_prep_fn;
        typedef std::function<void(const buf_t* src, uint_t nitem)> load_fill_fn;

        struct DistListFormat {
            const ListFormat m_local;
            /**
             * total number of items in the dataset
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

        protected:
            /**
             * this is kept protected so that only its two subclasses can be instantiated
             */
            DistListFormat(ItemFormat item_format, uint_t nitem_local,
                           uint_t nitem, uint_t nitem_displ, str_t leading_dim_name="item");
        };

        /**
         * use when the transaction involves non-overlapping ranges of items being allocated to each rank
         * suitable for Load and Save
         */
        struct PartDistListFormat : DistListFormat {
            PartDistListFormat(hdf5::dataset::ItemFormat item_format, uint_t nitem, str_t leading_dim_name="item");
        };

        /**
         * use when the transaction involves all items being handled on every rank
         * suitable for Load only
         */
        struct FullDistListFormat : DistListFormat {
            FullDistListFormat(hdf5::dataset::ItemFormat item_format, uint_t nitem, str_t leading_dim_name="item");
        };
    }
}

#endif //M7_DATASET_FORMAT_H