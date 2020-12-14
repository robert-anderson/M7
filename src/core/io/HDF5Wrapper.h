//
// Created by rja on 13/12/2020.
//

#ifndef M7_HDF5WRAPPER_H
#define M7_HDF5WRAPPER_H

#include <string>
#include <memory>
#include "src/core/parallel/MPIWrapper.h"
#include "src/core/nd/NdFormat.h"
#include "hdf5.h"

namespace hdf5 {

    static const std::array<hid_t, 10> types =
            {H5T_NATIVE_CHAR, H5T_NATIVE_SHORT, H5T_NATIVE_INT32, H5T_NATIVE_LONG,
             H5T_NATIVE_UCHAR, H5T_NATIVE_USHORT, H5T_NATIVE_UINT32, H5T_NATIVE_ULONG,
             H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE};

    template<typename T>
    static constexpr size_t type_ind() { return ~0ul; }

    template<> constexpr size_t type_ind<char>() { return 0;}
    template<> constexpr size_t type_ind<short int>() { return 1;}
    template<> constexpr size_t type_ind<int>() { return 2; }
    template<> constexpr size_t type_ind<long int>() { return 3; }
    template<> constexpr size_t type_ind<unsigned char>() { return 4; }
    template<> constexpr size_t type_ind<unsigned short int>() { return 5; }
    template<> constexpr size_t type_ind<unsigned int>() { return 6; }
    template<> constexpr size_t type_ind<unsigned long int>() { return 7; }
    template<> constexpr size_t type_ind<float>() { return 8; }
    template<> constexpr size_t type_ind<double>() { return 9; }

    template<typename T>
    const hid_t& type(){return types[type_ind<T>()];}

    struct Group;

    struct File {
        bool m_writemode;
        hid_t m_handle;
        File(std::string name, bool writemode): m_writemode(writemode){
            if (writemode) {
                m_handle = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                if (m_handle<0) mpi::stop_all("HDF5 file could not be opened for writing. It may be locked by another program");
            }
            else {
                if(!H5Fis_hdf5(name.c_str())) mpi::stop_all("Specified file is not HDF5 format");
                m_handle = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            }
        }

        ~File() {
            auto status = H5Fclose(m_handle);
            if (status) mpi::stop_all("HDF5 Error on closing file");
        }

        Group subgroup(std::string name);
    };

    template<typename T, size_t ndim>
    struct Array {
        NdFormat<ndim> m_format;
        size_t m_chunk_size;
        Array(NdFormat<ndim> format, size_t chunk_size):
        m_format(format), m_chunk_size(chunk_size){}

    };

    struct Group {
        bool m_writemode;
        hid_t m_handle;

        Group(hid_t parent, std::string name, bool writemode): m_writemode(writemode) {
            if (writemode) m_handle = H5Gcreate(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            else m_handle = H5Gopen2(parent, name.c_str(), H5P_DEFAULT);
        }

        Group subgroup(std::string name){
            return Group(m_handle, name, m_writemode);
        }

        ~Group(){
            auto status = H5Gclose(m_handle);
            if (status) mpi::stop_all("HDF5 Error on closing group");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const T& item, std::string name) {
            if (!m_writemode) mpi::stop_all("File is open for reading - save not permitted.");
            hsize_t shape = 1;
            auto dspace_handle = H5Screate_simple(1, &shape, NULL);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite (dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (const void*)&item);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            if (status) mpi::stop_all("HDF5 Error on primitive type save");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const std::complex<T>& item, std::string name) {
            if (!m_writemode) mpi::stop_all("File is open for reading - save not permitted.");
            hsize_t shape = 2;
            auto dspace_handle = H5Screate_simple(1, &shape, NULL);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (const void*)&item);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            if (status) mpi::stop_all("HDF5 Error on complex type save");
        }

//        template<typename T, size_t ndim>
//        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
//        save(const Array<T, ndim>& item, std::string name) {
//            if (!m_writemode) mpi::stop_all("File is open for reading - save not permitted.");
//            hsize_t shape = 1;
//            H5::DataSpace dataspace(1, &shape);
//            auto dataset = m_group.createDataSet(name, type<T>(), dataspace);
//            dataset.write((const void*)&item, type<T>(), dataspace);
//        }



        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        load(T& item, std::string name) {
            if (m_writemode) mpi::stop_all("File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)&item);
            if (status) mpi::stop_all("HDF5 Error on complex type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        load(const std::complex<T>& item, std::string name) const {
            if (m_writemode) mpi::stop_all("File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)&item);
            if (status) mpi::stop_all("HDF5 Error on primitive type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, T>::type
        get(std::string name){
            T tmp;
            load(tmp, name);
            return tmp;
        }
    };

}

#endif //M7_HDF5WRAPPER_H
