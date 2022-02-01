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
#include "hdf5_hl.h"

namespace hdf5 {

    static std::vector<hsize_t> convert_dims(const defs::inds &item_dims) {
        std::vector<hsize_t> out;
        out.reserve(item_dims.size());
        for (auto &i: item_dims) out.push_back(i);
        return out;
    }

#ifdef H5_HAVE_PARALLEL
    constexpr bool have_parallel = true;
#else
    constexpr bool have_parallel = false;
#endif

    static_assert(have_parallel, "HDF5 must be compiled with parallel functionality");

    static const std::array<hid_t, 11> types =
            {H5T_NATIVE_CHAR, H5T_NATIVE_SHORT, H5T_NATIVE_INT32, H5T_NATIVE_LONG,
             H5T_NATIVE_UCHAR, H5T_NATIVE_USHORT, H5T_NATIVE_UINT32, H5T_NATIVE_ULONG,
             H5T_NATIVE_ULLONG, H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE};

    template<typename T=void>
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
    constexpr size_t type_ind<unsigned long long int>() { return 8; }

    template<>
    constexpr size_t type_ind<float>() { return 9; }

    template<>
    constexpr size_t type_ind<double>() { return 10; }

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
                REQUIRE_NE(m_handle, 0,
                                "HDF5 file could not be opened for writing. It may be locked by another program");
            } else {
                REQUIRE_TRUE(H5Fis_hdf5(name.c_str()), "Specified file is not HDF5 format");
                m_handle = H5Fopen(name.c_str(), H5F_ACC_RDONLY, plist_id);
            }
            H5Pclose(plist_id);
        }

        ~File() {
            auto status = H5Fclose(m_handle);
            REQUIRE_TRUE(!status, "HDF5 Error on closing file");
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

    struct AttributeWriterBase {
        const hid_t m_parent_handle, m_h5type;
        const defs::inds m_shape;
        const size_t m_nelement;
        hid_t m_memspace_handle;
        hid_t m_handle;

        AttributeWriterBase(hid_t parent_handle, std::string name, const defs::inds &shape, hid_t h5type);

        ~AttributeWriterBase();

    private:
        void write_bytes(const char *src);

    public:
        template<typename T>
        static void write(hid_t parent, std::string name, const std::vector<T> &src) {
            AttributeWriterBase(parent, name, {src.size()}, type<T>()).write_bytes(
                    reinterpret_cast<const char *>(src.data()));
        }

        static void write(hid_t parent, std::string name, const std::string &src);

        static void write(hid_t parent, std::string name, const std::vector<std::string> &src);
    };


    struct AttributeReaderBase {
        const hid_t m_parent_handle, m_h5type;
        const defs::inds m_shape;
        const size_t m_nelement;
        hid_t m_handle;

        AttributeReaderBase(hid_t parent_handle, std::string name, const defs::inds &shape, hid_t h5type) :
                m_parent_handle(parent_handle), m_h5type(h5type), m_shape(shape),
                m_nelement(nd_utils::nelement(shape)) {
            m_handle = H5Aopen_name(m_parent_handle, name.c_str());
        }

        ~AttributeReaderBase() {
            H5Aclose(m_handle);
        }
    };

    template<typename T>
    struct AttributeReader : AttributeReaderBase {
        AttributeReader(hid_t parent_handle, std::string name, const defs::inds &shape) :
                AttributeReaderBase(parent_handle, name, shape, type<T>()) {}

    private:
        void read(char *dst) {
            auto status = H5Aread(m_handle, m_h5type, dst);
            DEBUG_ASSERT_TRUE(!status, "HDF5 attribute read failed");
        }

    public:
        void read(std::vector<T> &dst) {
            MPI_REQUIRE(dst.size() == m_nelement, "Destination length is incompatible with specified shape");
            read(reinterpret_cast<char *>(dst.data()));
        }
    };


    struct FileBase {
        static void check_is_hdf5(const std::string &name);

        const hid_t m_handle;
    protected:
        FileBase(hid_t handle) : m_handle(handle) {}

        ~FileBase() {
            H5Fclose(m_handle);
        }
    };

    struct FileWriter : FileBase {
        FileWriter(std::string name);
    };

    struct FileReader : FileBase {
        FileReader(std::string name);
    };

    /**
     * parent class for HDF5 group I/O
     *
     * there are three different categories of object whose I/O is handled by this HDF5 wrapper:
     *  1. Non-distributed primitive types and STL container specializations of primitive types
     *  2. Non-distributed user-defined classes
     *  3. Distributed user-defined classes (MappedTable)
     *
     * more complicated combinations of these are handled in the Archivable class, but the above are the building blocks.
     *
     * Categories 1 and 2 are typically for small metadata and verification information, writing these types requires
     * the identification of a "definitive rank", defaulting to the root, whose value of the written data is taken to be
     * the correct value to commit to disk.
     *
     * The most expensive category is 3 since this includes wavefunctions and multidimensional quantities. In recognition
     * of this, the HDF5 files are opened in collective mode, and so even operations in categories 1 and 2 (which could
     * be done in independent mode) involve null write operations on the non-definitive MPI rank, and all ranks read in
     * the same value.
     *
     * save methods defined in the Writer subclass, and load methods defined in the Reader subclass are for category 1,
     * i.e. the following types:
     *  primitive
     *  pointer to primitive (shaped)
     *  vector of primitives
     *  complex of primitives
     *  pointer to complex of primitive (shaped)
     *  vector of complex of primitives
     * these are not user-defined classes and therefore we cannot define save and load methods on these objects, so
     * these functions are handled in the subclass definitions below
     *
     */
    struct GroupBase {
        const std::string m_name;
        const hid_t m_parent_handle;
        const hid_t m_handle;
    protected:
        GroupBase(std::string name, hid_t parent_handle, hid_t handle);

        ~GroupBase();
    };

    /**
     * carries out all creation of datasets and Groups
     */
    struct GroupWriter : GroupBase {
        GroupWriter(std::string name, const FileWriter &parent);

        GroupWriter(std::string name, const GroupWriter &parent) :
                GroupBase(name, parent.m_handle,
                          H5Gcreate(parent.m_handle, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {}

        template<typename T>
        void write_attr(std::string name, const T &obj) {
            AttributeWriterBase::write(m_handle, name, obj);
        }

        /**
         * commit a single value of a primitive type (HDF5 scalar dataset) to disk
         * @tparam T
         *  primitive type (any type for which type_ind<T>() is not ~0ul)
         * @param name
         *  key in the HDF5 Group in which the value is to be stored
         * @param v
         *  value to store
         * @param irank
         *  index of MPI rank which stores the definitive value of v (other ranks write to null dataspace)
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(std::string name, const T &v, size_t irank=0ul) {
            auto dspace_handle = H5Screate(H5S_SCALAR);
            /**
             * make a null selection if this is not the rank we want to output the value of
             */
            if (!mpi::i_am(irank)) H5Sselect_none(dspace_handle);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, type<T>(), dspace_handle, dspace_handle, H5P_DEFAULT, static_cast<const void *>(&v));
            REQUIRE_FALSE_ALL(status, "HDF5 Error on primitive type save");
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
        }

        /**
         * commit a single value of a complex type (HDF5 simple dataset) to disk by writing a vector of size 2
         * @tparam T
         *  primitive type of the components of the complex number
         * @param name
         *  key in the HDF5 Group in which the value is to be stored
         * @param v
         *  value to store
         * @param irank
         *  index of MPI rank which stores the definitive value of v
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(std::string name, const std::complex<T> &v, size_t irank=0ul) {
            defs::inds shape = {2};
            save(name, reinterpret_cast<const T*>(&v), shape, {"real_imag"}, irank);
        }

        /**
         * commit a multidimensional array (HDF5 simple dataset) of a primitive type to disk
         * @tparam T
         *  primitive type of the elements of the array
         * @param name
         *  key in the HDF5 Group in which the value is to be stored
         * @param v
         *  pointer to the beginning of the array data
         * @param shape
         *  vector of dimensional extents
         * @param irank
         *  index of MPI rank which stores the definitive value of v
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(std::string name, const T* v, const defs::inds& shape, std::vector<std::string> dim_labels={}, size_t irank=0ul) {
            auto dims = convert_dims(shape);
            auto dspace_handle = H5Screate_simple(dims.size(), dims.data(), nullptr);
            /**
             * make a null selection if this is not the rank we want to output the value of
             */
            if (!mpi::i_am(irank)) H5Sselect_none(dspace_handle);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, type<T>(), dspace_handle, dspace_handle, H5P_DEFAULT, static_cast<const void*>(v));
            REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
            if (!dim_labels.empty()) {
                DEBUG_ASSERT_EQ(dim_labels.size(), dims.size(),
                                "Number of dim labels does not match number of dims");
                for (size_t idim = 0ul; idim < dims.size(); ++idim) {
                    H5DSset_label(dset_handle, idim, dim_labels[idim].c_str());
                    REQUIRE_FALSE(status, "HDF5 Error on dimension label assignment");
                }
            }
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
        }

        /**
         * save a multidimensional array of a complex type (HDF5 simple dataset) to disk by adding a new minor index of
         * extent 2
         * @tparam T
         *  primitive type of the complex components of elements of the array
         * @param name
         *  key in the HDF5 Group in which the value is to be stored
         * @param v
         *  pointer to the beginning of the array data
         * @param shape
         *  vector of dimensional extents
         * @param irank
         *  index of MPI rank which stores the definitive value of v
         * @return
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(std::string name, const std::complex<T>* v, const defs::inds& shape,
             std::vector<std::string> dim_labels={}, size_t irank=0ul) {
            dim_labels.push_back("real_imag");
            auto dims = shape;
            dims.push_back(2ul);
            save(name, reinterpret_cast<const T*>(v), dims, dim_labels, irank);
        }

        /**
         * wrapper for save of simple type in the case that the data source is a vector, but also has a multidimensional
         * shape
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(std::string name, const std::vector<T>& v, const defs::inds& shape, std::vector<std::string> dim_labels={}, size_t irank=0ul){
            REQUIRE_EQ_ALL(v.size(), nd_utils::nelement(shape), "vector and shape are incompatible");
            save(name, v.data(), shape, {}, irank);
        }

        /**
         * wrapper for save of vector in the case that it is not multidimensional
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(std::string name, const std::vector<T>& v, size_t irank=0ul){
            save(name, v, {v.size()}, {}, irank);
        }


        /**
         * save a vector of strings by creating a char array datatype for the longest element of the given vector v
         * @param name
         *  key in the HDF5 Group in which the value is to be stored
         * @param v
         *  vector of strings to store on disk
         * @param irank
         *  index of MPI rank which stores the definitive value of v
         */
        void save(std::string name, const std::vector<std::string>& v, size_t irank=0ul);

        /**
         * wrapper for save in the case that only a single string is to be stored
         */
        void save(std::string name, const std::string& v, size_t irank=0ul);
    };

    struct GroupReader : GroupBase {
        GroupReader(std::string name, const FileReader &parent);

        GroupReader(std::string name, const GroupReader &parent);

        bool child_exists(const std::string& name) const;

        size_t first_existing_child(const std::vector<std::string>& names) const;

        size_t nchild() const;

        std::string child_name(size_t ichild) const;

        int child_type(size_t i) const;

        std::vector<std::string> child_names(int type=-1) const;

    private:

        size_t get_dataset_ndim(std::string name);

        defs::inds get_dataset_shape(std::string name);

    public:

        /**
         * load a single value of a primitive type (HDF5 scalar dataset) from disk
         * @tparam T
         *  primitive type (any type for which type_ind<T>() is not ~0ul)
         * @param name
         *  key in the HDF5 Group in which the value is stored
         * @param v
         *  value to retrieve
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(std::string name, T& v) {
            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            auto dspace_handle = H5Screate(H5S_SCALAR);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
            auto status = H5Dread(dset_handle, type<T>(), dspace_handle, dspace_handle, H5P_DEFAULT, static_cast<void *>(&v));
            REQUIRE_FALSE_ALL(status, "HDF5 Error on primitive type load");
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
        }

        /**
         * load a single value of a complex type (HDF5 simple dataset) to disk by reinterpreting the complex number as
         * a length-2 array and reading a simple dataset from the file
         * @tparam T
         *  primitive type of the components of the complex number
         * @param name
         *  key in the HDF5 Group in which the value is stored
         * @param v
         *  value to store
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(std::string name, std::complex<T> &v, size_t irank=0ul) {
            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            load(name, reinterpret_cast<std::array<T, 2>&>(v)[0], irank);
            load(name, reinterpret_cast<std::array<T, 2>&>(v)[1], irank);
        }

        /**
         * load a multidimensional array (HDF5 simple dataset) of a primitive type from disk
         * @tparam T
         *  primitive type of the elements of the array
         * @param name
         *  key in the HDF5 Group in which the value is stored
         * @param v
         *  pointer to the beginning of the array data
         * @param shape
         *  vector of expected dimensional extents - throw error if this is not the same as the stored shape
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(std::string name, T* v, const defs::inds& shape){
            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            auto file_shape = get_dataset_shape(name);
            REQUIRE_EQ_ALL(shape, file_shape, "expected a container of a different shape");
            auto dims = convert_dims(shape);
            auto dspace_handle = H5Screate_simple(dims.size(), dims.data(), nullptr);
            auto dset_handle = H5Dopen1(m_handle, name.c_str());
            auto status = H5Dread(dset_handle, type<T>(), dspace_handle, dspace_handle, H5P_DEFAULT, static_cast<void*>(v));
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            REQUIRE_FALSE_ALL(status, "HDF5 Error on multidimensional load");
        }

        /**
         * load a multidimensional array (HDF5 simple dataset) of a complex type from disk
         * @tparam T
         *  primitive type of the real and imag components of elements of the array
         * @param name
         *  key in the HDF5 Group in which the value is stored
         * @param v
         *  pointer to the beginning of the array data
         * @param shape
         *  vector of expected dimensional extents - throw error if this is not the same as the stored shape
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(std::string name, std::complex<T>* v, const defs::inds& shape){
            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            auto complex_shape = shape;
            complex_shape.push_back(2ul);
            auto file_shape = get_dataset_shape(name);
            REQUIRE_EQ_ALL(complex_shape, file_shape, "expected a container of a different shape");
            auto dims = convert_dims(complex_shape);
            auto dspace_handle = H5Screate_simple(dims.size(), dims.data(), nullptr);
            auto dset_handle = H5Dopen1(m_handle, name.c_str());
            auto status = H5Dread(dset_handle, type<T>(), dspace_handle, dspace_handle, H5P_DEFAULT, static_cast<void*>(v));
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            REQUIRE_FALSE_ALL(status, "HDF5 Error on multidimensional load");
        }

        /**
         * convenient wrapper in the case that the destination is a vector but the source is shaped
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(std::string name, std::vector<T>& v, const defs::inds& shape){
            REQUIRE_EQ_ALL(v.size(), nd_utils::nelement(shape), "vector and shape are incompatible");
            load(name, v.data(), shape);
        }

        /**
         * convenient wrapper in the case that the destination is a vector
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(std::string name, std::vector<T>& v){
            load(name, v, {v.size()});
        }

        /**
         * convenient wrapper for scalar load
         */
        template<typename T>
        T load(std::string name) {
            T tmp;
            load(name, tmp);
            return tmp;
        }

        /**
         * convenient wrapper for vector load
         */
        template<typename T>
        std::vector<T> load_vector(std::string name) {
            auto nelement = nd_utils::nelement(get_dataset_shape(name));
            std::vector<T> tmp(nelement);
            load(name, tmp);
            return tmp;
        }

    };

    /**
     * Base class for the serialization of a distributed list of N-dimensional data.
     * Each rank handles (reads or writes) a number of items, and each item has a dimensionality which is constant
     * across all ranks. The list item dimension increases the overall dimensionality of the dataset by 1.
     */
    struct NdDistListBase {
        /**
         * HDF5 handle for the parent structure
         */
        const hid_t m_parent_handle;

        /**
         * shape of the list items
         */
        const std::vector<hsize_t> m_item_dims;
        /**
         * length of the item shape
         */
        const hsize_t m_ndim_item;
        /**
         * length of the list shape
         */
        const hsize_t m_ndim_list;
        /**
         * number of items the local rank is responsible for handling (reading / writing)
         */
        const hsize_t m_nitem_local;
        /**
         * total number of items handled over all ranks by this object
         */
        const hsize_t m_nitem_global;
        /**
         * largest number of items across all MPI ranks
         */
        const hsize_t m_nitem_local_max;
        /**
         * local shape of the list
         */
        const std::vector<hsize_t> m_list_dims_local;
        /**
         * global shape of the list
         */
        const std::vector<hsize_t> m_list_dims_global;
        /**
         * sum of m_nitem_local for all MPI ranks with index less than this rank
         */
        const hsize_t m_item_offset;

        /**
         * extent of the currently selected hyperslab in the HDF5 dataset
         */
        std::vector<hsize_t> m_hyperslab_counts;
        /**
         * multidimensional offset for the currently selected hyperslab
         */
        std::vector<hsize_t> m_hyperslab_offsets;
        /**
         * HDF5 handles
         */
        hid_t m_filespace_handle;
        hid_t m_dataset_handle;
        hid_t m_memspace_handle;
        /**
         * We use collective i/o mode, so when a rank has no reading or writing to do, we must point to a dummy memspace
         */
        hid_t m_none_memspace_handle;
        /**
         * datatype as an integer
         */
        hid_t m_h5type;
        /**
         * instance of object wrapping an HDF5 property list
         */
        CollectivePList m_coll_plist;

        /**
         * stick the local number of items onto the front of the item shape
         * @return
         *  list dims aka the overall shape of the local dataset
         */
        std::vector<hsize_t> get_list_dims_local();

        /**
         * stick the global number of items onto the front of the item shape
         * @return
         *  list dims aka the overall shape of the global dataset
         */
        std::vector<hsize_t> get_list_dims_global();

        /**
         * @return
         *  flat offset to the first item handled by this rank
         */
        hsize_t get_item_offset();

        NdDistListBase(hid_t parent_handle, std::string name, const defs::inds &item_dims, const size_t &nitem,
                       bool writemode, hid_t h5type);

        /**
         * sets the internal state to select the iitem-th item on this rank
         * @param iitem
         *  item index for selection
         */
        void select_hyperslab(const size_t &iitem);

        ~NdDistListBase();
    };

    /**
     * The writing case for a distributed N-dimensional list
     */
    struct NdDistListWriter : public NdDistListBase {
    private:
        NdDistListWriter(hid_t parent_handle, std::string name, const defs::inds &item_dims,
                         const size_t &nitem, hid_t h5type, const std::vector<std::string> &dim_labels = {});

    public:
        NdDistListWriter(FileWriter &parent, std::string name, const defs::inds &item_dims,
                         const size_t &nitem, hid_t h5type, const std::vector<std::string> &dim_labels = {}) :
                NdDistListWriter(parent.m_handle, name, item_dims, nitem, h5type, dim_labels) {}

        NdDistListWriter(GroupWriter &parent, std::string name, const defs::inds &item_dims,
                         const size_t &nitem, hid_t h5type, const std::vector<std::string> &dim_labels = {}) :
                NdDistListWriter(parent.m_handle, name, item_dims, nitem, h5type, dim_labels) {}

        /**
         * write the item stored at data to the iitem position within this rank
         * @param iitem
         *  item index determining the position on the HDF5 file buffer into which the data is copied
         * @param data
         *  pointer to data of any type. the counts in the hyperslab selection determine the number of elements n of the
         *  native type T corresponding to m_h5type are to be written. Thus the sizeof(T)*n bytes after data are copied
         */
        void write_h5item_bytes(const size_t &iitem, const void *data);

        template<typename T>
        void write_attr(std::string name, const T &obj) {
            AttributeWriterBase::write(m_dataset_handle, name, obj);
        }
    };

    /**
     * The reading case for a distributed N-dimensional list, splitting the reading load equally among all ranks, with
     * redistribution of items to be handled separately and subsquently
     */
    struct NdDistListReader : NdDistListBase {
        /**
         * interrogate the named dataset's object info for the number of shape elements
         * @param parent_handle
         *  HDF5 handle for group containing this dataset
         * @param name
         *  name of this dataset in group
         * @return
         *  length of the overall list shape
         */
        static size_t extract_list_ndim(hid_t parent_handle, std::string name) {
            auto status = H5Gget_objinfo(parent_handle, name.c_str(), 0, nullptr);
            REQUIRE_TRUE(!status, "Dataset \"" + name + "\" does not exist");
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            auto rank = H5Sget_simple_extent_ndims(dataspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            return rank;
        }
        /**
         * interrogate the named dataset's object info for the shape itself
         * @param parent_handle
         *  HDF5 handle for group containing this dataset
         * @param name
         *  name of this dataset in group
         * @return
         *  global list shape
         */
        static defs::inds extract_list_dims(hid_t parent_handle, std::string name) {
            auto rank = extract_list_ndim(parent_handle, name);
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

        /**
         * @param parent_handle
         *  HDF5 handle for group containing this dataset
         * @param name
         *  name of this dataset in group
         * @return
         *  shape of the list items themselves
         */
        static defs::inds extract_item_dims(hid_t parent_handle, std::string name) {
            auto list_dims = extract_list_dims(parent_handle, name);
            defs::inds out;
            out.reserve(list_dims.size() - 1);
            if (out.capacity()) out.insert(out.begin(), ++list_dims.cbegin(), list_dims.cend());
            return out;
        }

        /**
         * @param parent_handle
         *  HDF5 handle for group containing this dataset
         * @param name
         *  name of this dataset in group
         * @return
         *  just the global number of items
         */
        static size_t extract_nitem(hid_t parent_handle, std::string name) {
            return extract_list_dims(parent_handle, name)[0];
        }

        /**
         * share the elements out over the MPI ranks
         * @param parent_handle
         *  HDF5 handle for group containing this dataset
         * @param name
         *  name of this dataset in group
         * @return
         *  number of items this rank is responsible for reading int
         */
        static size_t local_nitem(hid_t parent_handle, std::string name) {
            auto nitem_tot = extract_nitem(parent_handle, name);
            auto nitem_share = nitem_tot / mpi::nrank();
            // give the remainder to the root rank
            if (mpi::i_am_root()) return nitem_share + nitem_tot % mpi::nrank();
            else return nitem_share;
        }

    private:
        NdDistListReader(hid_t parent_handle, std::string name, hid_t h5_type) :
                NdDistListBase(parent_handle, name, extract_item_dims(parent_handle, name),
                               local_nitem(parent_handle, name), false, h5_type) {}

    public:
        NdDistListReader(FileReader &parent, std::string name, hid_t h5_type) :
                NdDistListReader(parent.m_handle, name, h5_type) {}

        NdDistListReader(GroupReader &parent, std::string name, hid_t h5_type) :
                NdDistListReader(parent.m_handle, name, h5_type) {}

        /**
         * read the item stored at the iitem position into the data pointer
         * @param iitem
         *  item index determining the position on the HDF5 file buffer from which the data is copied
         * @param data
         *  pointer to data of any type. the counts in the hyperslab selection determine the number of elements n of the
         *  native type T corresponding to m_h5type are to be read.
         */
        void read_h5item_bytes(const size_t &iitem, void *data) {
            select_hyperslab(iitem);
            log::debug_("reading data...");
            if (data) {
                auto status = H5Dread(m_dataset_handle, m_h5type, m_memspace_handle,
                                      m_filespace_handle, m_coll_plist, data);
                DEBUG_ONLY(status);
                DEBUG_ASSERT_TRUE(!status, "HDF5 read failed");
            }
            else {
                auto status = H5Dread(m_dataset_handle, m_h5type, m_none_memspace_handle,
                                      m_filespace_handle, m_coll_plist, data);
                DEBUG_ONLY(status);
                DEBUG_ASSERT_TRUE(!status, "HDF5 read failed");
            }
            log::debug_("data read");
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
            REQUIRE_TRUE(!status, "HDF5 Error on closing group");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(const T &item, std::string name) {
            REQUIRE_TRUE(m_writemode, "File is open for reading - save not permitted.");
            hsize_t shape = 1;
            auto dspace_handle = H5Screate_simple(1, &shape, nullptr);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
            auto item_ptr = reinterpret_cast<const void*>(&item);
            auto status = H5Dwrite(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, item_ptr);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            REQUIRE_FALSE(status, "HDF5 Error on primitive type save");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(const std::complex<T> &item, std::string name) {
            REQUIRE_TRUE(m_writemode, "File is open for reading - save not permitted.");
            hsize_t shape = 2;
            auto dspace_handle = H5Screate_simple(1, &shape, nullptr);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), type<T>(), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
            auto item_ptr = reinterpret_cast<const void*>(&item);
            auto status = H5Dwrite(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, item_ptr);
            DEBUG_ONLY(status);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            REQUIRE_FALSE_ALL(status, "HDF5 Error on complex type save");
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
            REQUIRE_TRUE(!m_writemode, "File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto item_ptr = reinterpret_cast<void*>(&item);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, item_ptr);
            MPI_REQUIRE_ALL(!status, "HDF5 Error on complex type load");
        }

        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(const std::complex<T> &item, std::string name) const {
            REQUIRE_TRUE(!m_writemode, "File is open for writing - load not permitted.");
            auto dset_handle = H5Dopen2(m_handle, name.c_str(), H5P_DEFAULT);
            auto item_ptr = reinterpret_cast<void*>(&item);
            auto status = H5Dread(dset_handle, type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, item_ptr);
            REQUIRE_FALSE(status, "HDF5 Error on primitive type load");
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
