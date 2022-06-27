//
// Created by anderson on 27/06/2022.
//

#include "PropertyList.h"

hdf5::PList::PList(hid_t handle) : m_handle(handle){}

hdf5::PList::~PList() {
    H5Pclose(m_handle);
}

hdf5::PList::operator hid_t() const {
    return m_handle;
}

hdf5::AccessPList::AccessPList() : PList(H5Pcreate(H5P_FILE_ACCESS)) {
    H5Pset_fapl_mpio(m_handle, MPI_COMM_WORLD, MPI_INFO_NULL);
}

hdf5::CollectivePList::CollectivePList() : PList(H5Pcreate(H5P_DATASET_XFER)) {
    H5Pset_dxpl_mpio(m_handle, H5FD_MPIO_COLLECTIVE);
}
