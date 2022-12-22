//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_FILE_H
#define M7_HDF5_FILE_H

#include "PropertyList.h"
#include "NodeReader.h"
#include "NodeWriter.h"

namespace hdf5 {
    struct FileBase {
        const str_t m_fname;

        static bool is_hdf5(const str_t &fname);

        static void require_is_hdf5(const str_t &fname);

        FileBase(const str_t &fname) : m_fname(fname) {}
    };

    struct FileReader : NodeReader, FileBase {
    private:
        static hid_t get_handle(const str_t &fname) {
            require_is_hdf5(fname);
            AccessPList p_list;
            H5Pset_fapl_mpio(p_list.m_handle, MPI_COMM_WORLD, MPI_INFO_NULL);
            REQUIRE_TRUE(H5Fis_hdf5(fname.c_str()), "Specified file is not HDF5 format");
            return H5Fopen(fname.c_str(), H5F_ACC_RDONLY, p_list.m_handle);
        }

    public:
        FileReader(const str_t &fname) : NodeReader(get_handle(fname)), FileBase(fname) {}

        ~FileReader() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    };

    struct FileWriter : NodeWriter, FileBase {
    private:
        static hid_t get_handle(const str_t &fname) {
            AccessPList p_list;
            H5Pset_fapl_mpio(p_list, MPI_COMM_WORLD, MPI_INFO_NULL);
            auto handle = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, p_list);
            REQUIRE_NE(handle, 0, "HDF5 file could not be opened for writing. It may be locked by another program");
            return handle;
        }

    public:
        FileWriter(const str_t &fname) : NodeWriter(get_handle(fname)), FileBase(fname) {}

        ~FileWriter() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    };
}

#endif //M7_HDF5_FILE_H
