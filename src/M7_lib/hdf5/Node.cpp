//
// Created by anderson on 27/06/2022.
//

#include "Node.h"

hdf5::Node::Node(hid_t id) : m_id(id){}

hdf5::Node::operator hid_t() const {
    return m_id;
}

bool hdf5::Node::child_exists(const str_t& name) const {
    return H5Oexists_by_name(m_id, name.c_str(), H5P_DEFAULT) > 0;
}

H5O_info_t hdf5::get_object_info(hid_t obj_id) {
    H5O_info_t info;
    H5Oget_info(obj_id, &info);
    return info;
}
