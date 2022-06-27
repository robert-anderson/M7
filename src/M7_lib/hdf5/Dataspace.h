//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_DATASPACE_H
#define M7_HDF5_DATASPACE_H

#include "Type.h"

namespace hdf5 {

    struct DataSpace {
        const hid_t m_handle;
        const std::vector <hsize_t> m_shape;
        const hsize_t m_nelement;
    private:
        std::vector<hsize_t> make_shape() const;

    public:
        DataSpace(hid_t handle);

        DataSpace(const std::vector <hsize_t> &shape, bool select_none = false);

        ~DataSpace();

        operator hid_t() const;

        void select_none() const;
    };
}


#endif //M7_HDF5_DATASPACE_H
