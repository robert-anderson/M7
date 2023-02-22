//
// Created by anderson on 27/06/2022.
//

#include "Node.h"

hdf5::Node::Node(hid_t id) : m_id(id){}

hdf5::Node::operator hid_t() const {
    return m_id;
}

H5O_info_t hdf5::get_object_info(hid_t obj_id) {
    H5O_info_t info;
    H5Oget_info(obj_id, &info);
    return info;
}
