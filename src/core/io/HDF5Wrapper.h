//
// Created by rja on 13/12/2020.
//

#ifndef M7_HDF5WRAPPER_H
#define M7_HDF5WRAPPER_H

#include <string>
#include <memory>
#include "src/core/parallel/MPIWrapper.h"
#include "src/core/parallel/MPIAssert.h"
#include "src/core/nd/NdFormat.h"
#include "hdf5.h"

namespace hdf5 {

#ifdef H5_HAVE_PARALLEL
    constexpr bool have_parallel = true;
#else
    constexpr bool have_parallel = false;
#endif

    static_assert(have_parallel, "HDF5 must be compiled with parallel functionality");

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
            auto plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
            if (writemode) {
                m_handle = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                MPI_REQUIRE_ALL(m_handle>=0, "HDF5 file could not be opened for writing. It may be locked by another program");
            }
            else {
                MPI_REQUIRE_ALL(H5Fis_hdf5(name.c_str()), "Specified file is not HDF5 format");
                m_handle = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            }
            H5Pclose(plist_id);
        }

        ~File() {
            auto status = H5Fclose(m_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on closing file");
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
            MPI_REQUIRE_ALL(!status, "HDF5 Error on closing group");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const T& item, std::string name) {
            MPI_REQUIRE_ALL(m_writemode, "File is open for reading - save not permitted.");
            hsize_t shape = 1;
            auto dspace_handle = H5Screate_simple(1, &shape, NULL);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite (dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (const void*)&item);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on primitive type save");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        save(const std::complex<T>& item, std::string name) {
            MPI_REQUIRE_ALL(m_writemode, "File is open for reading - save not permitted.");
            hsize_t shape = 2;
            auto dspace_handle = H5Screate_simple(1, &shape, NULL);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (const void*)&item);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on complex type save");
        }

//        template<typename T, size_t ndim>
//        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
//        save(const Array<T, ndim>& item, std::string name) {
//            if (!m_writemode) mpi::abort("File is open for reading - save not permitted.");
//            hsize_t shape = 1;
//            H5::DataSpace dataspace(1, &shape);
//            auto dataset = m_group.createDataSet(name, type<T>(), dataspace);
//            dataset.write((const void*)&item, type<T>(), dataspace);
//        }



        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        load(T& item, std::string name) {
            MPI_REQUIRE_ALL(!m_writemode, "File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)&item);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on complex type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, void>::type
        load(const std::complex<T>& item, std::string name) const {
            MPI_REQUIRE_ALL(!m_writemode, "File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)&item);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on primitive type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>()!=~0ul, T>::type
        get(std::string name){
            T tmp;
            load(tmp, name);
            return tmp;
        }
    };

    template<typename T>
    struct VectorWriter {

        const hid_t m_parent;
        const std::string m_name;
        const hsize_t m_nblock;
        const hsize_t m_block_size;
        const hsize_t m_nblock_chunk;
        const hsize_t m_size;
        const hsize_t m_chunk_size;
        const hsize_t m_nchunk;
        hid_t m_plist;
        hid_t m_dataspace;
        hid_t m_dataset;
        hid_t m_chunk_space;

        size_t m_iblock = 0;
        size_t m_ichunk = 0;


        bool m_done = false;

        /**
         * move the chunk space window to the next position
         */
        void next_chunk(){
            const hsize_t ichunk_global = (m_ichunk++)*mpi::nrank()+mpi::irank();
            if (ichunk_global>=m_nchunk) m_done = true;
            if (m_done) return;

            const hsize_t start = ichunk_global*m_nblock_chunk;
            const hsize_t stride = 1;
            const hsize_t count = (ichunk_global<m_nchunk) ? m_nblock_chunk : (m_nblock - m_nblock_chunk*ichunk_global);
            ASSERT(count<=m_nblock_chunk);
            auto status = H5Sselect_hyperslab(m_chunk_space, H5S_SELECT_SET, &start, &stride, &count, &m_block_size);
            MPI_REQUIRE_ALL(status>=0, "HDF5 Error: could not select next hyperslab");
        }

        VectorWriter(hid_t parent, std::string name, size_t nblock, size_t block_size, size_t nblock_chunk):
        m_parent(parent), m_name(name), m_nblock(nblock), m_block_size(block_size),
        m_nblock_chunk(nblock_chunk), m_size(m_block_size*m_nblock),
        m_chunk_size(m_block_size*m_nblock_chunk), m_nchunk(integer_utils::divceil(m_size, m_chunk_size)){

            m_dataspace = H5Screate_simple(1, &m_chunk_size, nullptr);

            m_plist = H5Pcreate(H5P_DATASET_XFER);
            auto status = H5Pset_dxpl_mpio(m_plist, H5FD_MPIO_INDEPENDENT);
            MPI_REQUIRE_ALL(status>=0, "HDF5 Error: could not set data transfer mode to independent");

            m_dataset = H5Dcreate(parent, name.c_str(), type<T>(), m_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            m_chunk_space = H5Screate_simple(1, &m_chunk_size, nullptr);
            next_chunk();
        }

        ~VectorWriter() {
            H5Sclose(m_chunk_space);
            H5Dclose(m_dataset);
            H5Sclose(m_dataspace);
        }

        void write(const T* src, size_t nblock_write){
            auto status = H5Dwrite(m_dataset, type<T>(), m_chunk_space, m_dataspace, m_plist, src);
            MPI_REQUIRE_ALL(status<0, "HDF5 Error: vector write failed");
        }

    };

    template<typename T>
    class ListWriter {
        /*
         * parent HDF5 element (Group or File)
         */
        const hid_t m_parent;
        /*
         * dimensionality of an "item" e.g. if we are storing a list of matrices, this is 2
         */
        const hsize_t m_item_ndim;
        /*
         * shape of items
         */
        const std::vector<hsize_t> m_item_shape;
        /*
         * dimensionality of the list structure, e.g. if we are storing a list of matrices, this is 3
         */
        const hsize_t m_list_ndim;
        /*
         * 0th element is local number of rows
         */
        const std::vector<hsize_t> m_list_shape;
        /*
         * 0th element is total number of rows across all MPI ranks
         */
        const std::vector<hsize_t> m_list_shape_global;
        /*
         * each process has a hyperslab defining the location within the HDF5 file
         * counts: the number of T-values to be written / read
         * offsets: the position of the hyperslab (only the 0th value changes, rest are 0)
         */
        const std::vector<hsize_t> m_item_hyperslab_counts;
        std::vector<hsize_t> m_item_hyperslab_offsets;

        hid_t m_collective_write_plist_id;
        hid_t m_dataset_id;
        hid_t m_item_dataspace_id;
        hid_t m_filespace_id;

        std::vector<hsize_t> make_item_shape(const defs::inds& item_shape) {
            std::vector<hsize_t> out;
            out.reserve(item_shape.size());
            for (const auto &x: item_shape) out.push_back(x);
            return out;
        }

        std::vector<hsize_t> make_list_shape() {
            auto out = m_item_shape;
            out.insert(out.begin(), m_list_ndim);
            return out;
        }

        std::vector<hsize_t> make_list_shape_global() {
            auto out = m_item_shape;
            out.insert(out.begin(), mpi::all_sum(m_list_ndim));
            return out;
        }

        std::vector<hsize_t> make_item_hyperslab_counts() {
            std::vector<hsize_t> out;
            out.reserve(m_list_ndim);
            out.push_back(1);
            out.insert(++out.cbegin(), m_item_shape.cbegin(), m_item_shape.cend());
            return out;
        }

    public:

        ListWriter(hid_t parent, std::string name, defs::inds item_shape):
                m_parent(parent),
                m_item_ndim(item_shape.size()), m_item_shape(make_item_shape(item_shape)),
                m_list_ndim(m_item_ndim+1), m_list_shape(make_list_shape()),
                m_list_shape_global(make_list_shape_global()),
                m_item_hyperslab_counts(make_item_hyperslab_counts()),
                m_item_hyperslab_offsets(m_list_ndim, 0ul){
            /*
             * Create the dataspace for the dataset.
             */
            auto filespace = H5Screate_simple(m_list_ndim, m_list_shape_global.data(), NULL);
            /*
             * Create the dataset with default properties and close filespace.
             */
            m_dataset_id = H5Dcreate(m_parent, name.c_str(), type<T>(), filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            utils::print(m_item_hyperslab_counts);

            m_item_dataspace_id = H5Screate_simple(m_list_ndim, m_item_hyperslab_counts.data(), NULL);
            m_filespace_id = H5Dget_space(m_dataset_id);

            m_collective_write_plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(m_collective_write_plist_id, H5FD_MPIO_COLLECTIVE);
        }

        //ListWriter(const ListWriter<T>& other): ListWriter(other.m_parent){}

        ~ListWriter(){
            H5Dclose(m_dataset_id);
            H5Sclose(m_filespace_id);
            H5Sclose(m_item_dataspace_id);
            H5Pclose(m_collective_write_plist_id);
        }

        void write_item(const T* item){
            auto status = H5Dwrite(m_dataset_id, type<T>(), m_item_dataspace_id,
                              m_filespace_id, m_collective_write_plist_id, item);
            ASSERT(!status)
            ++m_item_hyperslab_offsets[0];
        }

    };

}

#endif //M7_HDF5WRAPPER_H
