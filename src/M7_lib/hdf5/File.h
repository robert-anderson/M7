//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_FILE_H
#define M7_HDF5_FILE_H

#include "PropertyList.h"
#include "Node.h"

namespace hdf5 {
    struct FileBase {
        const std::string m_fname;

        static bool is_hdf5(const std::string &fname);

        static void require_is_hdf5(const std::string &fname);

        FileBase(const std::string &fname) : m_fname(fname) {}
    };

    struct FileReader : NodeReader, FileBase {
    private:
        static hid_t get_handle(const std::string &fname) {
            require_is_hdf5(fname);
            AccessPList p_list;
            H5Pset_fapl_mpio(p_list.m_handle, MPI_COMM_WORLD, MPI_INFO_NULL);
            REQUIRE_TRUE(H5Fis_hdf5(fname.c_str()), "Specified file is not HDF5 format");
            return H5Fopen(fname.c_str(), H5F_ACC_RDONLY, p_list.m_handle);
        }

    public:
        FileReader(const std::string &fname) : NodeReader(get_handle(fname)), FileBase(fname) {}

        ~FileReader() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    };

    struct FileWriter : NodeWriter, FileBase {
    private:
        static hid_t get_handle(const std::string &fname) {
            AccessPList p_list;
            H5Pset_fapl_mpio(p_list, MPI_COMM_WORLD, MPI_INFO_NULL);
            auto handle = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, p_list);
            REQUIRE_NE(handle, 0, "HDF5 file could not be opened for writing. It may be locked by another program");
            return handle;
        }

    public:
        FileWriter(const std::string &fname) : NodeWriter(get_handle(fname)), FileBase(fname) {}

        ~FileWriter() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    };
}

#endif //M7_HDF5_FILE_H
