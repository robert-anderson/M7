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

    template<>
    constexpr size_t type_ind<char>() { return 0; }

    template<>
    constexpr size_t type_ind<short int>() { return 1; }

    template<>
    constexpr size_t type_ind<int>() { return 2; }

    template<>
    constexpr size_t type_ind<long int>() { return 3; }

    template<>
    constexpr size_t type_ind<unsigned char>() { return 4; }

    template<>
    constexpr size_t type_ind<unsigned short int>() { return 5; }

    template<>
    constexpr size_t type_ind<unsigned int>() { return 6; }

    template<>
    constexpr size_t type_ind<unsigned long int>() { return 7; }

    template<>
    constexpr size_t type_ind<float>() { return 8; }

    template<>
    constexpr size_t type_ind<double>() { return 9; }

    template<typename T>
    const hid_t &type() { return types[type_ind<T>()]; }

    struct Group;

    struct File {
        bool m_writemode;
        hid_t m_handle;

        File(std::string name, bool writemode) : m_writemode(writemode) {
            auto plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

            if (writemode) {
                m_handle = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                MPI_REQUIRE_ALL(m_handle >= 0,
                                "HDF5 file could not be opened for writing. It may be locked by another program");
            } else {
                MPI_REQUIRE_ALL(H5Fis_hdf5(name.c_str()), "Specified file is not HDF5 format");
                m_handle = H5Fopen(name.c_str(), H5F_ACC_RDONLY, plist_id);
            }
            H5Pclose(plist_id);
        }

        ~File() {
            auto status = H5Fclose(m_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on closing file");
        }

        Group subgroup(std::string name);
    };

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

    struct FileBase {
        static void check_is_hdf5(const std::string &name) {
            MPI_REQUIRE(H5Fis_hdf5(name.c_str()), "Specified file is not HDF5 format");
        }

        const hid_t m_handle;
    protected:
        FileBase(hid_t handle) : m_handle(handle) {}

        ~FileBase() {
            H5Fclose(m_handle);
        }
    };

    struct FileWriter : FileBase {
        FileWriter(std::string name) : FileBase(H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, AccessPList())) {
            MPI_REQUIRE_ALL(m_handle >= 0,
                            "HDF5 file could not be opened for writing. It may be locked by another program");
        }
    };

    struct FileReader : FileBase {
        FileReader(std::string name) : FileBase(
                (check_is_hdf5(name), H5Fopen(name.c_str(), H5F_ACC_RDONLY, AccessPList()))) {
            MPI_REQUIRE_ALL(m_handle >= 0, "HDF5 file could not be opened for reading.");
        }
    };

    struct GroupBase {
        const hid_t m_parent_handle;
        const hid_t m_handle;
    protected:
        GroupBase(hid_t parent_handle, hid_t handle) :
                m_parent_handle(parent_handle), m_handle(handle) {}

        ~GroupBase() {
            auto status = H5Gclose(m_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on closing group");
        }
    };

    struct GroupWriter : GroupBase {
        GroupWriter(std::string name, const FileWriter &parent) :
                GroupBase(parent.m_handle,
                          H5Gcreate(parent.m_handle, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {}

        GroupWriter(std::string name, const GroupWriter &parent) :
                GroupBase(parent.m_handle,
                          H5Gcreate(parent.m_handle, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {}
    };

    struct GroupReader : GroupBase {
        GroupReader(std::string name, const FileReader &parent) :
                GroupBase(parent.m_handle, H5Gopen2(parent.m_handle, name.c_str(), H5P_DEFAULT)) {}

        GroupReader(std::string name, const GroupReader &parent) :
                GroupBase(parent.m_handle, H5Gopen2(parent.m_handle, name.c_str(), H5P_DEFAULT)) {}
    };


    struct NdListBase {
        const hid_t m_parent_handle;

        const std::vector<hsize_t> m_item_dims;
        const hsize_t m_ndim_item;
        const hsize_t m_ndim_list;
        const hsize_t m_nitem_local;
        const hsize_t m_nitem_global;
        const std::vector<hsize_t> m_list_dims_local;
        const std::vector<hsize_t> m_list_dims_global;
        const hsize_t m_item_offset;

        std::vector<hsize_t> m_hyperslab_counts;
        std::vector<hsize_t> m_hyperslab_offsets;
        hid_t m_filespace_handle;
        hid_t m_dataset_handle;
        hid_t m_memspace_handle;
        hid_t m_h5type;
        CollectivePList m_coll_plist;

        std::vector<hsize_t> get_item_dims(const defs::inds &item_dims) {
            std::vector<hsize_t> out;
            out.reserve(item_dims.size());
            for (auto &i: item_dims) out.push_back(i);
            return out;
        }

        std::vector<hsize_t> get_list_dims_local() {
            std::vector<hsize_t> out;
            out.reserve(m_ndim_list);
            out.push_back(m_nitem_local);
            out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
            return out;
        }

        std::vector<hsize_t> get_list_dims_global() {
            std::vector<hsize_t> out;
            out.reserve(m_ndim_list);
            out.push_back(m_nitem_global);
            out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
            return out;
        }

        hsize_t get_item_offset() {
            std::vector<hsize_t> tmp(mpi::nrank());
            mpi::all_gather(m_nitem_local, tmp);
            hsize_t out = 0ul;
            for (size_t irank = 0ul; irank < mpi::irank(); ++irank) out += tmp[irank];
            return out;
        }

        NdListBase(hid_t parent_handle, std::string name, const defs::inds &item_dims, const size_t &nitem,
                   bool writemode, hid_t h5type) :
                m_parent_handle(parent_handle),
                m_item_dims(get_item_dims(item_dims)),
                m_ndim_item(item_dims.size()),
                m_ndim_list(item_dims.size() + 1),
                m_nitem_local(nitem),
                m_nitem_global(mpi::all_sum(m_nitem_local)),
                m_list_dims_local(get_list_dims_local()),
                m_list_dims_global(get_list_dims_global()),
                m_item_offset(get_item_offset()),
                m_hyperslab_counts(m_ndim_list, 0ul),
                m_hyperslab_offsets(m_ndim_list, 0ul),
                m_h5type(h5type){
            m_filespace_handle = H5Screate_simple(m_ndim_list, m_list_dims_global.data(), nullptr);

            /*
             * Create the dataset with default properties and close filespace.
             */
            if (writemode)
                m_dataset_handle = H5Dcreate(m_parent_handle, name.c_str(), m_h5type, m_filespace_handle,
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            else
                m_dataset_handle = H5Dopen1(m_parent_handle, name.c_str());

            H5Sclose(m_filespace_handle);

            m_filespace_handle = H5Dget_space(m_dataset_handle);

            m_hyperslab_counts = m_list_dims_local;
            // select one item at a time
            m_hyperslab_counts[0] = 1;
            m_memspace_handle = H5Screate_simple(m_ndim_list, m_hyperslab_counts.data(), nullptr);
        }

        void select_hyperslab(const size_t& iitem){
            MPI_ASSERT(iitem < m_nitem_local, "iitem exceeds local limit");
            m_hyperslab_offsets[0] = m_item_offset + iitem;
            H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                nullptr, m_hyperslab_counts.data(), nullptr);
        }

        ~NdListBase() {
            H5Sclose(m_filespace_handle);
            H5Dclose(m_dataset_handle);
            H5Sclose(m_memspace_handle);
        }

    };


    struct NdListWriterBase : public NdListBase {
    private:
        NdListWriterBase(hid_t parent_handle, std::string name, const defs::inds &item_dims, const size_t &nitem, hid_t h5type) :
                NdListBase(parent_handle, name, item_dims, nitem, true, h5type) {}

    public:
        NdListWriterBase(FileWriter &parent, std::string name, const defs::inds &item_dims, const size_t &nitem, hid_t h5type) :
                NdListWriterBase(parent.m_handle, name, item_dims, nitem, h5type) {}

        NdListWriterBase(GroupWriter &parent, std::string name, const defs::inds &item_dims, const size_t &nitem, hid_t h5type) :
                NdListWriterBase(parent.m_handle, name, item_dims, nitem, h5type) {}

        void write_h5item_bytes(const size_t& iitem, const void *data) {
            select_hyperslab(iitem);
            auto status = H5Dwrite(m_dataset_handle, m_h5type, m_memspace_handle,
                                   m_filespace_handle, m_coll_plist, data);
            MPI_ASSERT(!status, "HDF5 write failed");
        }
    };

    template<typename T>
    struct NdListWriter : public NdListWriterBase {
        NdListWriter(FileWriter &parent, std::string name, const defs::inds &item_dims, const size_t &nitem) :
                NdListWriterBase(parent.m_handle, name, item_dims, nitem, type<T>()) {}

        NdListWriter(GroupWriter &parent, std::string name, const defs::inds &item_dims, const size_t &nitem) :
                NdListWriterBase(parent.m_handle, name, item_dims, nitem, type<T>()) {}

        void write_h5item(const size_t& iitem, const T *data) {
            write_h5item_bytes(iitem, data);
        }
    };


    struct NdListReaderBase : NdListBase {

        static size_t extract_list_rank(hid_t parent_handle, std::string name) {
            auto status = H5Gget_objinfo(parent_handle, name.c_str(), 0, nullptr);
            MPI_REQUIRE(!status, "Dataset \"" + name + "\" does not exist");
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            auto rank = H5Sget_simple_extent_ndims(dataspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return rank;
        }

        static defs::inds extract_list_dims(hid_t parent_handle, std::string name) {
            auto rank = extract_list_rank(parent_handle, name);
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            std::vector<hsize_t> dims(rank, 0ul);
            H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            defs::inds out;
            out.reserve(dims.size());
            for (const auto &i: dims) out.push_back(i);
            return out;
        }

        static defs::inds extract_item_dims(hid_t parent_handle, std::string name) {
            auto list_dims = extract_list_dims(parent_handle, name);
            defs::inds out;
            out.reserve(list_dims.size() - 1);
            if (out.capacity()) out.insert(out.begin(), ++list_dims.cbegin(), list_dims.cend());
            return out;
        }

        static size_t extract_nitem(hid_t parent_handle, std::string name) {
            return extract_list_dims(parent_handle, name)[0];
        }

        /*
         * share the elements out over the MPI ranks
         */
        static size_t rank_nitem(hid_t parent_handle, std::string name) {
            auto nitem_tot = extract_nitem(parent_handle, name);
            auto nitem_share = nitem_tot/mpi::nrank();
            // give the remainder to the root rank
            if (mpi::i_am_root()) return nitem_share + nitem_tot%mpi::nrank();
            else return nitem_share;
        }

    private:
        NdListReaderBase(hid_t parent_handle, std::string name, hid_t h5_type) :
                NdListBase(parent_handle, name, extract_item_dims(parent_handle, name),
                           rank_nitem(parent_handle, name), false, h5_type) {}

    public:
        NdListReaderBase(FileReader &parent, std::string name, hid_t h5_type) :
                NdListReaderBase(parent.m_handle, name, h5_type) {}

        NdListReaderBase(GroupReader &parent, std::string name, hid_t h5_type) :
                NdListReaderBase(parent.m_handle, name, h5_type) {}

        void read_h5item_bytes(const size_t& iitem, void *data) {
            select_hyperslab(iitem);
            auto status = H5Dread(m_dataset_handle, m_h5type, m_memspace_handle,
                                  m_filespace_handle, m_coll_plist, data);
            MPI_ASSERT(!status, "HDF5 read failed");
        }
    };


    template<typename T>
    struct NdListReader : public NdListReaderBase {
        NdListReader(FileWriter &parent, std::string name) :
                NdListReaderBase(parent.m_handle, name, type<T>()) {}

        NdListReader(GroupWriter &parent, std::string name) :
                NdListReaderBase(parent.m_handle, name, type<T>()) {}

        void read_h5item(const size_t& iitem, T *data) {
            read_h5item_bytes(iitem, data);
        }
    };


    template<typename T, size_t ndim>
    struct Array {
        NdFormat<ndim> m_format;
        size_t m_chunk_size;

        Array(NdFormat<ndim> format, size_t chunk_size) :
                m_format(format), m_chunk_size(chunk_size) {}

    };


    struct Group {
        bool m_writemode;
        hid_t m_handle;

        Group(hid_t parent, std::string name, bool writemode) : m_writemode(writemode) {
            if (writemode) m_handle = H5Gcreate(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            else m_handle = H5Gopen2(parent, name.c_str(), H5P_DEFAULT);
        }

        Group subgroup(std::string name) {
            return Group(m_handle, name, m_writemode);
        }

        ~Group() {
            auto status = H5Gclose(m_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on closing group");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(const T &item, std::string name) {
            MPI_REQUIRE_ALL(m_writemode, "File is open for reading - save not permitted.");
            hsize_t shape = 1;
            auto dspace_handle = H5Screate_simple(1, &shape, nullptr);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (const void *) &item);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on primitive type save");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(const std::complex<T> &item, std::string name) {
            MPI_REQUIRE_ALL(m_writemode, "File is open for reading - save not permitted.");
            hsize_t shape = 2;
            auto dspace_handle = H5Screate_simple(1, &shape, nullptr);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (const void *) &item);
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
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(T &item, std::string name) {
            MPI_REQUIRE_ALL(!m_writemode, "File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *) &item);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on complex type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(const std::complex<T> &item, std::string name) const {
            MPI_REQUIRE_ALL(!m_writemode, "File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *) &item);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on primitive type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, T>::type
        get(std::string name) {
            T tmp;
            load(tmp, name);
            return tmp;
        }
    };


}

#endif //M7_HDF5WRAPPER_H
