//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_GROUP_H
#define M7_HDF5_GROUP_H

#include "NodeReader.h"
#include "NodeWriter.h"

namespace hdf5 {

    struct GroupWriter : NodeWriter {
        GroupWriter(const NodeWriter &node, str_t name);

        ~GroupWriter();
    };

    struct GroupReader : NodeReader {
        GroupReader(const NodeReader &node, str_t name);

        ~GroupReader();
    };
}

#endif //M7_HDF5_GROUP_H
