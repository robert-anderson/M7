//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_DATASET_H
#define M7_HDF5_DATASET_H

#include "Dataspace.h"

namespace hdf5 {

    struct DatasetReader {
        const DataSpace m_space;
        const hid_t m_handle;

        static uint_t get_ndim(hid_t parent_handle, const std::string& name);

        template<typename T>
        static std::vector<T> get_shape(hid_t parent_handle, const std::string& name) {
            auto ndim = get_ndim(parent_handle, name);
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            REQUIRE_GT_ALL(dataset, 0, log::format("no such dataset \"{}\"", name));
            auto dataspace = H5Dget_space(dataset);
            std::vector<hsize_t> shape(ndim);
            H5Sget_simple_extent_dims(dataspace, shape.data(), nullptr);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return convert::vector<T>(shape);
        }

        static uint_t get_nelement(hid_t parent_handle, const std::string& name);

    public:
        const Type m_type;

        DatasetReader(hid_t parent_handle, const std::string& name);

        ~DatasetReader();

        void read(void* dst) const;
    };

    struct DatasetWriter {
        const DataSpace m_space;
        const hid_t m_handle;
        const std::vector<std::string> m_dim_names;
    public:
        const Type m_type;

        DatasetWriter(hid_t parent_handle, const std::string& name, const std::vector<hsize_t>& shape, Type type,
                      std::vector<std::string> dim_names={}, uint_t irank=0ul);

        ~DatasetWriter();

        void write(const void* src) const;
    };
}


#endif //M7_HDF5_DATASET_H
