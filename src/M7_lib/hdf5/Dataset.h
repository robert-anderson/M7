//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_DATASET_H
#define M7_HDF5_DATASET_H

#include "Dataspace.h"

namespace hdf5 {

    struct DatasetReader {
        const DataSpace m_space;
        const hid_t m_handle;

        static uint_t get_ndim(hid_t parent_handle, const str_t& name);

        static v_t<hsize_t> get_hdf5_shape(hid_t parent_handle, const str_t& name);

        static uint_t get_nelement(hid_t parent_handle, const str_t& name);

    public:
        const Type m_type;

        DatasetReader(hid_t parent_handle, const str_t& name);

        ~DatasetReader();

        void read(void* dst) const;
    };

    struct DatasetWriter {
        const DataSpace m_space;
        const hid_t m_handle;
        const strv_t m_dim_names;
    public:
        const Type m_type;

        DatasetWriter(hid_t parent_handle, const str_t& name, const v_t<hsize_t>& shape, Type type,
                      strv_t dim_names={}, uint_t irank=0ul);

        ~DatasetWriter();

        void write(const void* src) const;
    };
}


#endif //M7_HDF5_DATASET_H
