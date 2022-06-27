//
// Created by anderson on 27/06/2022.
//

#include "Dataset.h"

uint_t hdf5::DatasetReader::get_ndim(const hdf5::NodeReader &node, const std::string &name) {
    auto status = H5Gget_objinfo(node.m_handle, name.c_str(), 0, nullptr);
    REQUIRE_TRUE(!status, "Dataset \"" + name + "\" does not exist");
    auto dataset = H5Dopen1(node.m_handle, name.c_str());
    auto dataspace = H5Dget_space(dataset);
    auto rank = H5Sget_simple_extent_ndims(dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return rank;
}

uint_t hdf5::DatasetReader::get_nelement(const hdf5::NodeReader &node, const std::string &name) {
    return nd::nelement(get_shape<hsize_t>(node, name));
}

hdf5::DatasetReader::DatasetReader(const hdf5::NodeReader &node, const std::string &name) :
        m_space(get_shape<hsize_t>(node, name)), m_handle(H5Dopen1(node.m_handle, name.c_str())),
        m_type(H5Dget_type(m_handle)) {}

hdf5::DatasetReader::~DatasetReader() {
    H5Dclose(m_handle);
}

void hdf5::DatasetReader::read(char *dst) const {
    auto status = H5Dread(m_handle, m_type, m_space, m_space, H5P_DEFAULT, dst);
    REQUIRE_FALSE_ALL(status, "HDF5 Error on dataset load");
}

hdf5::DatasetWriter::DatasetWriter(const hdf5::NodeWriter &node, const std::string &name,
                                   const std::vector<hsize_t> &shape, Type type, std::vector<std::string> dim_names,
                                   uint_t irank) :
        m_space(DataSpace(shape, !mpi::i_am(irank))),
        m_handle(H5Dcreate2(node.m_handle, name.c_str(), type, m_space.m_handle, H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT)), m_dim_names(std::move(dim_names)), m_type(H5Dget_type(m_handle)) {}

hdf5::DatasetWriter::~DatasetWriter() {
    H5Dclose(m_handle);
}

void hdf5::DatasetWriter::write(const char *src) const {
    auto status = H5Dwrite(m_handle, m_type, m_space, m_space, H5P_DEFAULT, src);
    REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
    if (!m_dim_names.empty()) {
        DEBUG_ASSERT_EQ(m_dim_names.size(), m_space.m_shape.size(),
                        "Number of dim labels does not match number of dims");
        for (uint_t idim = 0ul; idim < m_space.m_shape.size(); ++idim) {
            H5DSset_label(m_handle, idim, m_dim_names[idim].c_str());
            REQUIRE_FALSE(status, "HDF5 Error on dimension label assignment");
        }
    }
}
