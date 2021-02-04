//
// Created by rja on 13/12/2020.
//

#include "HDF5Wrapper.h"

hdf5::Group hdf5::File::subgroup(std::string name) {
    return Group(m_handle, name, m_writemode);
}
