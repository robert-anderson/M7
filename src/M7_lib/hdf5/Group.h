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
        GroupWriter(const NodeWriter &node, std::string name);

        ~GroupWriter();
    };

    struct GroupReader : NodeReader {
        GroupReader(const NodeReader &node, std::string name);

        ~GroupReader();
    };
}

#endif //M7_HDF5_GROUP_H
