//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_FILE_H
#define M7_HDF5_FILE_H

#include "Node.h"

namespace hdf5 {
    struct FileBase {
        const str_t m_fname;

        static bool is_hdf5(const str_t &fname);

        static void require_is_hdf5(const str_t &fname);

        FileBase(const str_t &fname) : m_fname(fname) {}
    };

    struct FileReader : NodeReader, FileBase {
    private:
        static hid_t get_handle(const str_t &fname);

    public:
        FileReader(const str_t &fname) : NodeReader(get_handle(fname)), FileBase(fname) {}

        ~FileReader() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    };

    struct FileWriter : NodeWriter, FileBase {
    private:
        static hid_t get_handle(const str_t &fname);

    public:
        FileWriter(const str_t &fname) : NodeWriter(get_handle(fname)), FileBase(fname) {}

        ~FileWriter() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
        }
    };
}

#endif //M7_HDF5_FILE_H
