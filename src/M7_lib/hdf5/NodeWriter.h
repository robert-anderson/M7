//
// Created by rja on 22/12/22.
//

#ifndef M7_NODEWRITER_H
#define M7_NODEWRITER_H

#include "M7_lib/util/Pointer.h"
#include "Node.h"

namespace hdf5 {

    struct NodeWriter : Node {
        NodeWriter(hid_t handle): Node(handle){}

        void save_attr(const Attr& attr) const;

        template<typename T>
        void save_attr(const str_t& name, const T& v) const {
            save_attr({v, name});
        }
    };
}


#endif //M7_NODEWRITER_H
