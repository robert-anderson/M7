//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_DATASET_H
#define M7_HDF5_DATASET_H

#include "Node.h"

namespace hdf5 {

    struct DatasetReader {
        const DataSpace m_space;
        const hid_t m_handle;

        static uint_t get_ndim(const NodeReader& node, const std::string& name) {
            auto status = H5Gget_objinfo(node.m_handle, name.c_str(), 0, nullptr);
            REQUIRE_TRUE(!status, "Dataset \"" + name + "\" does not exist");
            auto dataset = H5Dopen1(node.m_handle, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            auto rank = H5Sget_simple_extent_ndims(dataspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return rank;
        }

        template<typename T>
        static std::vector<T> get_shape(const NodeReader& node, const std::string& name) {
            auto ndim = get_ndim(node.m_handle, name);
            auto dataset = H5Dopen1(node.m_handle, name.c_str());
            REQUIRE_GT_ALL(dataset, 0, log::format("no such dataset \"{}\"", name));
            auto dataspace = H5Dget_space(dataset);
            std::vector<hsize_t> shape(ndim);
            H5Sget_simple_extent_dims(dataspace, shape.data(), nullptr);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return convert::vector<T>(shape);
        }

        static uint_t get_nelement(const NodeReader& node, const std::string& name) {
            return nd::nelement(get_shape<hsize_t>(node, name));
        }

    public:
        const hid_t m_h5type;

        DatasetReader(const NodeReader& node, const std::string& name) :
                m_space(get_shape<hsize_t>(node, name)),
                m_handle(H5Dopen1(node.m_handle, name.c_str())),
                m_h5type(H5Dget_type(m_handle)) {}

        ~DatasetReader() {
            H5Dclose(m_handle);
        }

        void read_bytes(char* dst) const {
            auto status = H5Dread(m_handle, m_h5type, m_space, m_space, H5P_DEFAULT, dst);
            REQUIRE_FALSE_ALL(status, "HDF5 Error on dataset load");
        }

        template<typename T>
        void read(T* dst) {
            REQUIRE_TRUE(types_equal<T>(m_h5type), "element type is at odds with the stored type");
            read_bytes(reinterpret_cast<char*>(dst));
        }

        template<typename T>
        void read(std::complex<T>* dst) {
            REQUIRE_TRUE(types_equal<T>(m_h5type), "element type is at odds with the stored type");
            REQUIRE_EQ(m_space.m_shape.back(), 2ul, "complex arrays should have a minor extent of 2");
            read_bytes(reinterpret_cast<char*>(dst));
        }
    };

    struct DatasetWriter {
        const DataSpace m_space;
        const hid_t m_handle;
        const std::vector<std::string> m_dim_names;
    public:
        const hid_t m_h5type;

        DatasetWriter(const NodeWriter& node, const std::string& name,
                      const std::vector<hsize_t>& shape, hid_t h5type,
                      std::vector<std::string> dim_names={}, uint_t irank=0ul):
                m_space(DataSpace(shape, !mpi::i_am(irank))),
                m_handle(H5Dcreate2(node.m_handle, name.c_str(), h5type, m_space.m_handle, H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT)), m_dim_names(std::move(dim_names)), m_h5type(H5Dget_type(m_handle)) {}

        ~DatasetWriter(){
            H5Dclose(m_handle);
        }

        void write_bytes(const char* src) const {
            auto status = H5Dwrite(m_handle, m_h5type, m_space, m_space, H5P_DEFAULT, src);
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

        template<typename T>
        void write(const T* src) const {
            REQUIRE_TRUE(types_equal<T>(m_h5type), "element type is at odds with the stored type");
            write_bytes(reinterpret_cast<const char*>(src));
        }

        template<typename T>
        void write(const std::complex<T>* src) {
            REQUIRE_TRUE(types_equal<T>(m_h5type), "element type is at odds with the stored type");
            REQUIRE_EQ(m_space.m_shape.back(), 2ul, "complex arrays should have a minor extent of 2");
            write_bytes(reinterpret_cast<char*>(src));
        }
    };
}


#endif //M7_HDF5_DATASET_H
