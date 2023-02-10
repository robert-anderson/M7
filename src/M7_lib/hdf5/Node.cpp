//
// Created by anderson on 27/06/2022.
//

#include "Node.h"

hdf5::Node::Node(hid_t handle) : m_handle(handle){}

hdf5::Node::operator hid_t() const {
    return m_handle;
}

H5O_info_t hdf5::get_object_info(hid_t obj_handle) {
    H5O_info_t info;
    H5Oget_info(obj_handle, &info);
    return info;
}
