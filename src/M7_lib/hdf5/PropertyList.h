//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_PROPERTYLIST_H
#define M7_HDF5_PROPERTYLIST_H

#include "Type.h"

namespace hdf5 {
    struct PList {
        const hid_t m_handle;

        PList(hid_t handle);

        ~PList();

        operator hid_t() const;
    };

    struct AccessPList : PList {
        AccessPList();
    };

    struct CollectivePList : PList {
        CollectivePList();
    };
}


#endif //M7_HDF5_PROPERTYLIST_H
