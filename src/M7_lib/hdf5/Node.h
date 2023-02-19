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

        bool child_exists(const str_t& name) const {
            return H5Oexists_by_name(m_handle, name.c_str(), H5P_DEFAULT);
        }
    };

    H5O_info_t get_object_info(hid_t obj_handle);

    struct NodeReader : Node {
        NodeReader(hid_t handle) : Node(handle) {}

        template<typename T>
        T load_attr(const str_t& name) const {
            T v;
            auto success = Attr(m_handle, name).parse(v);
            REQUIRE_TRUE(success, "HDF5 attribute load failed without default value");
            return v;
        }

        template<typename T>
        T load_attr(const str_t& name, T default_) const {
            T v;
            Attr(m_handle, name).parse(v, default_);
            return v;
        }
    };

    struct NodeWriter : Node {
        NodeWriter(hid_t handle): Node(handle){}

        void save_attr(const Attr& attr) const {
            attr.save(m_handle);
        }

        template<typename T>
        void save_attr(const str_t& name, const T& v) const {
            save_attr({v, name});
        }
    };


    struct GroupWriter : NodeWriter {
        GroupWriter(const NodeWriter &node, str_t name) :
            NodeWriter(H5Gcreate(node, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {}

        ~GroupWriter() {
            auto status = H5Gclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing group");
        }
    };

    struct GroupReader : NodeReader {
        GroupReader(const NodeReader &node, str_t name) :
            NodeReader(H5Gopen2(node, name.c_str(), H5P_DEFAULT)) {}

        ~GroupReader() {
            auto status = H5Gclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing group");
        }
    };
}

#endif //M7_HDF5_NODE_H