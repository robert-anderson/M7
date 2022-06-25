//
// Created by anderson on 24/06/2022.
//

#include "test_core/defs.h"
#include "M7_lib/io/HDF5Wrapper.h"

namespace hdf5_new {

    static std::vector<hsize_t> convert_dims(const uintv_t &item_dims) {
        std::vector<hsize_t> out;
        out.reserve(item_dims.size());
        for (auto &i: item_dims) out.push_back(i);
        return out;
    }

#ifdef H5_HAVE_PARALLEL
    constexpr bool c_have_parallel = true;
#else
    constexpr bool c_have_parallel = false;
#endif

    static_assert(c_have_parallel, "HDF5 must be compiled with parallel functionality");

    static const std::array<hid_t, 12> types =
            {0, H5T_NATIVE_CHAR, H5T_NATIVE_SHORT, H5T_NATIVE_INT32, H5T_NATIVE_LONG,
                    H5T_NATIVE_UCHAR, H5T_NATIVE_USHORT, H5T_NATIVE_UINT32, H5T_NATIVE_ULONG,
                    H5T_NATIVE_ULLONG, H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE};

    template<typename T=void>
    static constexpr uint_t type_ind() { return 0; }

    template<>
    constexpr uint_t type_ind<char>() { return 1; }

    template<>
    constexpr uint_t type_ind<short int>() { return 2; }

    template<>
    constexpr uint_t type_ind<int>() { return 3; }

    template<>
    constexpr uint_t type_ind<long int>() { return 4; }

    template<>
    constexpr uint_t type_ind<unsigned char>() { return 5; }

    template<>
    constexpr uint_t type_ind<unsigned short int>() { return 6; }

    template<>
    constexpr uint_t type_ind<unsigned int>() { return 7; }

    template<>
    constexpr uint_t type_ind<unsigned long int>() { return 8; }

    template<>
    constexpr uint_t type_ind<unsigned long long int>() { return 9; }

    template<>
    constexpr uint_t type_ind<float>() { return 10; }

    template<>
    constexpr uint_t type_ind<double>() { return 11; }

    template<typename T>
    const hid_t &type() {
        typedef arith::comp_t<T> comp_t;
        static_assert(type_ind<comp_t>(), "type has no HDF5 equivalent");
        return types[type_ind<comp_t>()];
    }

    hsize_t type_size(hid_t h5type) {
        return H5Tget_size(h5type);
    }

    struct PList {
        hid_t m_handle;

        ~PList() {
            H5Pclose(m_handle);
        }

        operator hid_t() const {
            return m_handle;
        }
    };

    struct AccessPList : PList {
        AccessPList() : PList{H5Pcreate(H5P_FILE_ACCESS)} {
            H5Pset_fapl_mpio(m_handle, MPI_COMM_WORLD, MPI_INFO_NULL);
        }
    };

    struct CollectivePList : PList {
        CollectivePList() : PList{H5Pcreate(H5P_DATASET_XFER)} {
            H5Pset_dxpl_mpio(m_handle, H5FD_MPIO_COLLECTIVE);
        }
    };


    struct Node {
        const hid_t m_handle;
        Node(hid_t handle): m_handle(handle){}
    };

    struct Attr {
        const hid_t m_handle;
        const hid_t m_h5type;
        const hsize_t m_size;
        const hsize_t m_nelement;
    private:
        std::vector<char> m_buf;
    protected:
        char* m_buf_ptr;
    public:
        Attr(hid_t handle, hid_t h5type):
            m_handle(handle), m_h5type(h5type), m_size(H5Aget_storage_size(m_handle)),
            m_nelement(m_size / type_size(m_h5type)), m_buf(m_size, 0), m_buf_ptr(&m_buf[0]){
            DEBUG_ASSERT_EQ(m_nelement* type_size(m_h5type), m_size,
                            "returned storage size is not an integral multiple of element size");
        }
        Attr(hid_t handle): Attr(handle, H5Aget_type(handle)){}

        ~Attr(){
            H5Aclose(m_handle);
        }

        static bool exists(hid_t handle, const std::string& name) {
            return H5Aexists(handle, name.c_str()) >= 0;
        }
    };

    struct AttrReader : Attr {
    private:
        static hid_t get_handle(const Node& node, const std::string& name) {
            REQUIRE_TRUE(exists(node.m_handle, name), "can't create AttrReader object: not found");
            return H5Aopen(node.m_handle, name.c_str(), AccessPList());
        }
    public:
        AttrReader(const Node& node, const std::string& name): Attr(get_handle(node, name)){}


        template<typename T>
        void read(T& v, uint_t ielement=0ul) const {
            REQUIRE_EQ(sizeof(T), type_size(m_h5type), "incompatible size for read");
            REQUIRE_LT(ielement, m_nelement, "element index OOB");
            H5Aread(m_handle, m_h5type, m_buf_ptr);
            v = reinterpret_cast<const T*>(m_buf_ptr)[ielement];
        }

        template<typename T>
        T read(uint_t ielement=0ul) const {
            T v;
            read(v, ielement);
            return v;
        }

        template<typename T>
        void read(std::vector<T>& v) const {
            v.clear();
            v.reserve(m_nelement);
            for (uint_t i=0; i<m_nelement; ++i) v.push_back(read<T>(i));
        }
    };

    struct ReadNode : Node {
    protected:
        H5O_info_t m_info;
    public:
        ReadNode(hid_t handle): Node(handle){
            H5Oget_info(m_handle, &m_info);
        }
    };

    struct WriteNode : Node {
        WriteNode(hid_t handle): Node(handle){}
    };




    struct ReadFile : ReadNode {
    private:
        static hid_t get_handle(const std::string& fname) {
            AccessPList p_list;
            H5Pset_fapl_mpio(p_list.m_handle, MPI_COMM_WORLD, MPI_INFO_NULL);
            REQUIRE_TRUE(H5Fis_hdf5(fname.c_str()), "Specified file is not HDF5 format");
            return H5Fopen(fname.c_str(), H5F_ACC_RDONLY, p_list.m_handle);
        }
    public:
        ReadFile(const std::string& fname): ReadNode(get_handle(fname)) {}

        ~ReadFile() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }

//                m_handle = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//                REQUIRE_NE(m_handle, 0,
//                           "HDF5 file could not be opened for writing. It may be locked by another program");
    };

    /*
    struct FileWriter : NodeWriter {
        FileWriter(std::string fname) {
            auto plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
            REQUIRE_TRUE(H5Fis_hdf5(fname.c_str()), "Specified file is not HDF5 format");
            m_handle = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_id);
            H5Pclose(plist_id);
        }

        FileReaderer() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    }
     */


//    struct NodeWriter : Node {
//        NodeWriter(hid_t handle): Node(handle){}
//    };

//    struct FileWriter : NodeWriter {
//
//    };

}

TEST(Hdf5FcidumpReader, Header) {
    using namespace hdf5_new;
    const std::string fname = PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5";
    ReadFile file(fname);
    AttrReader attr_reader(file, "ORBSYM");
    std::cout << attr_reader.m_size << std::endl;
    uintv_t v;
    attr_reader.read(v);
    std::cout << v << std::endl;



//    auto err = H5Aexists_by_name(file.m_handle, ".", "CORE_ENERGY", hdf5::AccessPList());
//    std::cout << err << std::endl;
}