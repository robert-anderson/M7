//
// Created by rja on 13/12/2020.
//

#ifndef M7_HDF5WRAPPER_H
#define M7_HDF5WRAPPER_H

#include <string>
#include <memory>
#include "hdf5/serial/H5Cpp.h"
#include "src/core/parallel/MPIWrapper.h"
#include "src/core/nd/NdFormat.h"

namespace hdf5 {

    static const std::array<H5::PredType, 10> types =
            {H5::PredType::NATIVE_CHAR, H5::PredType::NATIVE_SHORT, H5::PredType::NATIVE_INT32, H5::PredType::NATIVE_LONG,
             H5::PredType::NATIVE_UCHAR, H5::PredType::NATIVE_USHORT, H5::PredType::NATIVE_UINT32, H5::PredType::NATIVE_ULONG,
             H5::PredType::NATIVE_FLOAT, H5::PredType::NATIVE_DOUBLE};

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
    const H5::PredType& type(){return types[type_ind<T>()];}

    struct Group;

    struct File {
        bool m_writemode;
        std::unique_ptr<H5::H5File> m_file;
        File(std::string name, bool writemode);

        Group group(std::string name);
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
        H5::Group m_group;

        Group subgroup(std::string name){
            return {m_writemode, m_writemode ? m_group.createGroup(name) : m_group.openGroup(name)};
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const T& item, std::string name) {
            if (!m_writemode) mpi::stop_all("File is open for reading - save not permitted.");
            hsize_t shape = 1;
            H5::DataSpace dataspace(1, &shape);
            auto dataset = m_group.createDataSet(name, type<T>(), dataspace);
            dataset.write((const void*)&item, type<T>(), dataspace);
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const std::complex<T>& item, std::string name) {
            if (!m_writemode) mpi::stop_all("File is open for reading - save not permitted.");
            hsize_t shape = 2;
            H5::DataSpace dataspace(1, &shape);
            auto dataset = m_group.createDataSet(name, type<T>(), dataspace);
            dataset.write((const void*)&item, type<T>(), dataspace);
        }

        template<typename T, size_t ndim>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const Array<T, ndim>& item, std::string name) {
            if (!m_writemode) mpi::stop_all("File is open for reading - save not permitted.");
            hsize_t shape = 1;
            H5::DataSpace dataspace(1, &shape);
            auto dataset = m_group.createDataSet(name, type<T>(), dataspace);
            dataset.write((const void*)&item, type<T>(), dataspace);
        }


        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        load(T& item, std::string name) const {
            if (m_writemode) mpi::stop_all("File is open for writing - load not permitted.");
            hsize_t shape = 1;
            H5::DataSpace dataspace(1, &shape);
            auto dataset = m_group.openDataSet(name);
            dataset.read((void*)&item, type<T>(), dataspace);
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        load(std::complex<T>& item, std::string name) const {
            if (m_writemode) mpi::stop_all("File is open for writing - load not permitted.");
            hsize_t shape = 2;
            H5::DataSpace dataspace(1, &shape);
            auto dataset = m_group.openDataSet(name);
            dataset.read((void*)&item, type<T>(), dataspace);
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, T>::type
        get(std::string name){
            T tmp;
            load(tmp, name);
            return tmp;
        }
    };

    struct Array {

    };

}


#endif //M7_HDF5WRAPPER_H
