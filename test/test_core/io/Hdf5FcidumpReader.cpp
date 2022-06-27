//
// Created by anderson on 24/06/2022.
//

#include "test_core/defs.h"
#include "M7_lib/io/HDF5Wrapper.h"
#include "M7_lib/io/FcidumpTextFileReader.h"

using namespace hdf5;
namespace hdf5_dev {
    struct DatasetReader {
        const DataSpace m_space;
        const hid_t m_handle;

        uint_t get_dataset_ndim(hid_t parent_handle, const std::string& name) const {
            auto status = H5Gget_objinfo(parent_handle, name.c_str(), 0, nullptr);
            REQUIRE_TRUE(!status, "Dataset \"" + name + "\" does not exist");
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            auto rank = H5Sget_simple_extent_ndims(dataspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return rank;
        }

        uintv_t get_dataset_shape(hid_t parent_handle, const std::string& name) const {
            auto ndim = get_dataset_ndim(parent_handle, name);
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            REQUIRE_GT_ALL(dataset, 0, log::format("no such dataset \"{}\"", name));
            auto dataspace = H5Dget_space(dataset);
            std::vector<hsize_t> dims(ndim, 0ul);
            H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            uintv_t out;
            out.reserve(dims.size());
            for (const auto& i: dims) out.push_back(i);
            return out;
        }

    public:
        const hid_t m_h5type;

        //const hsize_t m_nelement;
        DatasetReader(hid_t parent_handle, const std::string& name) :
                m_space(convert::vector<hsize_t>(get_dataset_shape(parent_handle, name))),
                m_handle(H5Dopen1(parent_handle, name.c_str())),
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

        template<typename T>
        void read(std::string* dst) {
            //todo
//            REQUIRE_TRUE(types_equal<T>(m_h5type), "element type is at odds with the stored type");
//            read_bytes(reinterpret_cast<char*>(dst));
        }
    };
}

TEST(Hdf5FcidumpReader, Header) {
    hdf5::FileReader fr(PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5");
    auto s = fr.read_attr<std::string>("MOLCAS_MODULE", std::string("123"));
    std::cout << s << "|" << std::endl;

    hdf5_dev::DatasetReader dr(fr.m_handle, "FOCK_INDEX");
    std::cout << dr.m_space.m_shape << std::endl;

    std::vector<int64_t> inds(12, 0);
    dr.read(inds.data(), inds.size());

    std::cout << inds << std::endl;

//    FcidumpInfo info(fr);
}