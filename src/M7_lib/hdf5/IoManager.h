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

        struct ItemFormat {
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

            ItemFormat(Type h5_type, uintv_t shape, strv_t dim_names, bool add_complex_dim);

            ItemFormat() : ItemFormat(Type(), {}, {}, false){}
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

            ListFormat(ItemFormat item_format, uint_t nitem);
        };
        /**
         * the save-load functions only need access to the list format, not the details of the distribution
         */
        typedef std::function<const void*(uint_t i, const ListFormat& format, uint_t max_nitem_per_op)> save_fn;
        typedef std::function<void*(uint_t i, const ListFormat& format, uint_t max_nitem_per_op)> load_fn;

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
            DistListFormat(ItemFormat item_format, uint_t nitem_local, uint_t nitem, uint_t nitem_displ);
        };

        /**
         * use when the transaction involves non-overlapping ranges of items being allocated to each rank
         * suitable for Load and Save
         */
        struct PartDistListFormat : DistListFormat {
            PartDistListFormat(hdf5::dataset::ItemFormat item_format, uint_t nitem) :
                DistListFormat(item_format, nitem, mpi::all_sum(nitem),
                    mpi::counts_to_displs_consec(mpi::all_gathered(nitem))[mpi::irank()]){}
        };

        /**
         * use when the transaction involves all items being handled on every rank
         * suitable for Load only
         */
        struct FullDistListFormat : DistListFormat {
            FullDistListFormat(hdf5::dataset::ItemFormat item_format, uint_t nitem) :
                DistListFormat(item_format, nitem, nitem, 0ul){}
        };
    }
}

#endif //M7_HDF5_BUFFERER_H