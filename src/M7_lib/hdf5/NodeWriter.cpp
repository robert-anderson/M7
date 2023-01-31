//
// Created by rja on 22/12/22.
//

#include "M7_lib/util/Vector.h"
#include "NodeWriter.h"

void hdf5::NodeWriter::save_attr(const hdf5::Attr& attr) const {
    attr.save(m_handle);
}