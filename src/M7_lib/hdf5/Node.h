//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_NODE_H
#define M7_HDF5_NODE_H

#include "Attr.h"

namespace hdf5 {

    struct Node {
        const hid_t m_handle;
        Node(hid_t handle);
        operator hid_t() const;
    };

    static constexpr uint_t c_default_max_nitem_per_op = 16000000;
}

#endif //M7_HDF5_NODE_H