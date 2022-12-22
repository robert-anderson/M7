//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_NODE_H
#define M7_HDF5_NODE_H

#include "Attr.h"
#include "Dataset.h"
#include "IoManager.h"
#include "M7_lib/util/Pointer.h"

namespace hdf5 {

    struct Node {
        const hid_t m_handle;
        Node(hid_t handle);
        operator hid_t() const;
        bool attr_exists(const str_t& name) const;
    };
}

#endif //M7_HDF5_NODE_H