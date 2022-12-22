//
// Created by anderson on 27/06/2022.
//

#include "Node.h"

hdf5::Node::Node(hid_t handle) : m_handle(handle){}

hdf5::Node::operator hid_t() const {
    return m_handle;
}

bool hdf5::Node::attr_exists(const str_t& name) const {
    return H5Aexists(m_handle, name.c_str());
}