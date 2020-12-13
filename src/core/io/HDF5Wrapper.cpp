//
// Created by rja on 13/12/2020.
//

#include "HDF5Wrapper.h"

hdf5::File::File(std::string name, bool writemode) : m_writemode(writemode),
                                                     m_file(new H5::H5File(name, writemode ? H5F_ACC_TRUNC : H5F_ACC_RDONLY)){}

hdf5::Group hdf5::File::group(std::string name) {
    return {m_writemode, m_writemode ? m_file->createGroup(name) : m_file->openGroup(name)};
}
