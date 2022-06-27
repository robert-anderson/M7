//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_GROUP_H
#define M7_HDF5_GROUP_H

#include "Node.h"

namespace hdf5 {
    /**
     * carries out all creation of datasets and Groups
     */
    struct GroupWriter : NodeWriter {
        GroupWriter(const NodeWriter &node, std::string name) :
                NodeWriter(H5Gcreate(node, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {}
    };

    struct GroupReader : NodeReader {
        GroupReader(const NodeReader &node, std::string name) :
                NodeReader(H5Gopen2(node, name.c_str(), H5P_DEFAULT)) {}
    };
}

#endif //M7_HDF5_GROUP_H
