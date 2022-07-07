//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_NODE_H
#define M7_HDF5_NODE_H

#include "Attr.h"
#include "Dataset.h"

namespace hdf5 {

    struct Node {
        const hid_t m_handle;
        Node(hid_t handle);
        operator hid_t() const;
        bool attr_exists(const str_t& name) const;
    };


    struct NodeReader : Node {
    protected:
        H5O_info_t m_info;
    public:
        NodeReader(hid_t handle) : Node(handle) {
            H5Oget_info(m_handle, &m_info);
        }

    private:
        template<typename T>
        void read_attr_fn(const str_t& name, T& v, T default_) const {
            if (!attr_exists(name)) {
                v = default_;
                return;
            }
            AttrReader attr(m_handle, name);
            attr.read(&v);
        }

        template<typename T>
        void read_attr_fn(const str_t& name, v_t<T>& v, v_t<T> default_) const {
            if (!attr_exists(name)) {
                v = default_;
                return;
            }
            AttrReader attr(m_handle, name);
            v.resize(attr.m_nelement);
            attr.read(v.data());
        }
    public:

        uintv_t dataset_shape(const str_t& name) const {
            return convert::vector<uint_t>(DatasetReader::get_hdf5_shape(m_handle, name));
        }

        uint_t dataset_nelement(const str_t& name) const {
            return DatasetReader::get_nelement(m_handle, name);
        }

        template<typename T>
        void read_data(const str_t& name, T *v, uint_t size) const {
            DatasetReader dr(*this, name);
            REQUIRE_EQ(dr.m_type, Type(v), "incorrect type");
            REQUIRE_EQ(size, dr.m_space.m_nelement, "number of elements does not match read size");
            dr.read(v);
        }

        template<typename T>
        void read_data(const str_t& name, std::complex<T> *v, uint_t size) const {
            read_data(name, &arith::real_ref(*v), size*2);
        }

        template<typename T>
        void read_data(const str_t& name, T &v) const {
            read_data(name, &v, 1);
        }

        template<typename T>
        void read_data(const str_t& name, v_t<T> &v) const {
            auto nelement = DatasetReader::get_nelement(m_handle, name);
            v.resize(nelement);
            read_data(name, v.data(), v.size());
        }

        template<typename T>
        void read_data(const str_t& name, v_t<std::complex<T>> &v) const {
            auto nelement = DatasetReader::get_nelement(m_handle, name);
            REQUIRE_FALSE(nelement&1ul, "complex read requires even number of elements");
            v.resize(nelement/2);
            read_data(name, v.data(), v.size());
        }

        template<typename T>
        T read_data(const str_t& name) const {
            T v{};
            read_data(name, v);
            return v;
        }

        template<typename T>
        T read_attr(const str_t& name, T default_ = {}) const {
            T v;
            read_attr_fn(name, v, default_);
            return v;
        }

        bool child_exists(const str_t& name) const;

        uint_t first_existing_child(const strv_t& names) const;

        uint_t nchild() const;

        str_t child_name(uint_t ichild) const;

        int child_type(uint_t i) const;

        strv_t child_names(int type=-1) const;


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
        load(const str_t& name, T& v) const {
            REQUIRE_TRUE(0, "deprecated");
            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            auto dspace_handle = H5Screate(H5S_SCALAR);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), Type(v), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
            auto status = H5Dread(dset_handle, Type(v), dspace_handle, dspace_handle, H5P_DEFAULT, &v);
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
        load(const str_t& name, std::complex<T> &v) const {
            REQUIRE_TRUE(0, "deprecated");

            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            load(name, reinterpret_cast<std::array<T, 2>&>(v)[0]);
            load(name, reinterpret_cast<std::array<T, 2>&>(v)[1]);
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
        load(const str_t& name, T* v, const uintv_t& shape){
            REQUIRE_TRUE(0, "deprecated");

            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            auto file_shape = get_dataset_shape(name);
            REQUIRE_EQ_ALL(shape, file_shape, "expected a container of a different shape");
            auto dims = convert::vector<hsize_t>(shape);
            auto dspace_handle = H5Screate_simple(dims.size(), dims.data(), nullptr);
            auto dset_handle = H5Dopen1(m_handle, name.c_str());
            auto status = H5Dread(dset_handle, Type(v), dspace_handle, dspace_handle, H5P_DEFAULT, v);
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
        load(const str_t& name, std::complex<T>* v, const uintv_t& shape){
            REQUIRE_TRUE(0, "deprecated");

            DEBUG_ASSERT_TRUE(child_exists(name), "Can't read from non-existent object");
            auto complex_shape = shape;
            complex_shape.push_back(2ul);
            auto file_shape = get_dataset_shape(name);
            REQUIRE_EQ_ALL(complex_shape, file_shape, "expected a container of a different shape");
            auto dims = convert::vector<hsize_t>(shape);
            auto dspace_handle = H5Screate_simple(dims.size(), dims.data(), nullptr);
            auto dset_handle = H5Dopen1(m_handle, name.c_str());
            auto status = H5Dread(dset_handle, Type(v), dspace_handle, dspace_handle, H5P_DEFAULT, v);
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
            REQUIRE_FALSE_ALL(status, "HDF5 Error on multidimensional load");
        }

        /**
         * convenient wrapper in the case that the destination is a vector but the source is shaped
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(const str_t& name, v_t<T>& v, const uintv_t& shape){
            REQUIRE_TRUE(0, "deprecated");

            REQUIRE_EQ_ALL(v.size(), nd::nelement(shape), "vector and shape are incompatible");
            load(name, v.data(), shape);
        }

        /**
         * convenient wrapper in the case that the destination is a vector
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        load(const str_t& name, v_t<T>& v){
            load(name, v, {v.size()});
        }

        /**
         * convenient wrapper for scalar load
         */
        template<typename T>
        T load(const str_t& name) const {
            T tmp;
            load(name, tmp);
            return tmp;
        }

        /**
         * convenient wrapper for vector load
         */
        template<typename T>
        v_t<T> load_vector(const str_t& name) const {
            auto nelement = nd::nelement(get_dataset_shape(name));
            v_t<T> tmp(nelement);
            load(name, tmp);
            return tmp;
        }

    private:

        uint_t get_dataset_ndim(const str_t& name) const;

        uintv_t get_dataset_shape(const str_t& name) const;
    };

    struct NodeWriter : Node {
        NodeWriter(hid_t handle): Node(handle){}

        template<typename T>
        void write_attr(const str_t& name, const T& v) const {
            AttrWriter attr(m_handle, name, {1}, Type(&v));
            attr.write(&v, 1);
        }

        template<typename T>
        void write_attr(const str_t& name, const v_t<T>& v) const {
            AttrWriter attr(m_handle, name, {v.size()}, Type(v.data()));
            attr.write(v.data(), v.size());
        }

        void write_attr(const str_t& name, const str_t& v) const {
            AttrWriter attr(m_handle, name, {1}, Type(&v));
            attr.write(v.c_str(), 1);
        }

        void write_attr(const str_t& name, const strv_t& v) const {
            AttrWriter attr(m_handle, name, {v.size()}, Type(v.data()));
            attr.write(v.data()->c_str(), 1);
        }

        /**
         * since this is a public method, the shape type is uintv_t, which is converted to the vector type used internally within the hdf5 namespace
         * @tparam T
         * @param name
         * @param v
         * @param shape
         * @param dim_names
         */
        template<typename T>
        void write_data(const str_t& name, const T *v, const uintv_t& shape, const strv_t& dim_names={}, uint_t irank=0) {
            using namespace convert;
            DatasetWriter(*this, name, vector<hsize_t>(shape), Type(v), dim_names, irank).write(v);
        }
        template<typename T>
        void write_data(const str_t& name, const std::complex<T> *v, const uintv_t& shape, const strv_t& dim_names={}, uint_t irank=0) {
            auto tmp_shape = shape;
            tmp_shape.push_back(2);
            auto tmp_dim_names = dim_names;
            tmp_dim_names.push_back("real_imag");
            write_data(name, &arith::real_ref(*v), tmp_shape, tmp_dim_names, irank);
        }
        template<typename T>
        void write_data(const str_t& name, const v_t<T> &v, const uintv_t& shape, const strv_t& dim_names={}, uint_t irank=0) {
            write_data(name, v.data(), shape, dim_names, irank);
        }
        template<typename T>
        void write_data(const str_t& name, const T &v, uint_t irank=0) {
            write_data(name, &v, {1}, {}, irank);
        }
        template<typename T>
        void write_data(const str_t& name, const v_t<T> &v, uint_t irank=0) {
            write_data(name, v.data(), {v.size()}, {}, irank);
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
        save(const str_t& name, const T &v, uint_t irank=0ul) const {
            REQUIRE_TRUE(0, "deprecated");
            auto dspace_handle = H5Screate(H5S_SCALAR);
            /**
             * make a null selection if this is not the rank we want to output the value of
             */
            if (!mpi::i_am(irank)) H5Sselect_none(dspace_handle);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), Type(v), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, Type(v), dspace_handle, dspace_handle, H5P_DEFAULT, &v);
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
        save(const str_t& name, const std::complex<T> &v, uint_t irank=0ul) const {
            uintv_t shape = {2};
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
        save(const str_t& name, const T* v, const uintv_t& shape, strv_t dim_labels={}, uint_t irank=0ul) const {
            REQUIRE_TRUE(0, "deprecated");
            auto dims = convert::vector<hsize_t>(shape);
            auto dspace_handle = H5Screate_simple(dims.size(), dims.data(), nullptr);
            /**
             * make a null selection if this is not the rank we want to output the value of
             */
            if (!mpi::i_am(irank)) H5Sselect_none(dspace_handle);
            auto dset_handle = H5Dcreate2(m_handle, name.c_str(), Type(v), dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            auto status = H5Dwrite(dset_handle, Type(v), dspace_handle, dspace_handle, H5P_DEFAULT, v);
            REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
            if (!dim_labels.empty()) {
                DEBUG_ASSERT_EQ(dim_labels.size(), dims.size(),
                                "Number of dim labels does not match number of dims");
                for (uint_t idim = 0ul; idim < dims.size(); ++idim) {
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
        save(const str_t& name, const std::complex<T>* v, const uintv_t& shape,
             strv_t dim_labels={}, uint_t irank=0ul) const {

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
        save(const str_t& name, const v_t<T>& v, const uintv_t& shape, strv_t dim_labels={}, uint_t irank=0ul) const {
            REQUIRE_EQ_ALL(v.size(), nd::nelement(shape), "vector and shape are incompatible");
            save(name, v.data(), shape, dim_labels, irank);
        }

        /**
         * wrapper for save of vector in the case that it is not multidimensional
         */
        template<typename T>
        typename std::enable_if<type_ind<T>() != ~0ul, void>::type
        save(const str_t& name, const v_t<T>& v, uint_t irank=0ul) const {
            save(name, v, {v.size()}, {}, irank);
        }


        /**
         * save a vector of strings by creating a char array dtype for the longest element of the given vector v
         * @param name
         *  key in the HDF5 Group in which the value is to be stored
         * @param v
         *  vector of strings to store on disk
         * @param irank
         *  index of MPI rank which stores the definitive value of v
         */
        void save(const str_t& name, const strv_t& v, uint_t irank=0ul) const {
            REQUIRE_TRUE(0, "deprecated");
            auto memtype = H5Tcopy (H5T_C_S1);
            auto longest = std::max_element(
                    v.cbegin(), v.cend(),[](const str_t& s1, const str_t& s2){return s1.size()<s2.size();});
            auto size = longest->size();
            v_t<char> buffer(size*v.size());
            uint_t istr = 0ul;
            for (auto& str: v) std::strcpy(buffer.data()+(istr++)*size, str.c_str());

            auto status = H5Tset_size(memtype, size);
            REQUIRE_FALSE_ALL(status, "failed to create string type");

            /*
             * Create dataset with a null dataspace.
             */
            v_t<hsize_t> shape = {v.size()};
            auto dspace_handle = H5Screate_simple(1, shape.data(), NULL);
            /**
             * make a null selection if this is not the rank we want to output the value of
             */
            if (!mpi::i_am(irank)) H5Sselect_none(dspace_handle);
            auto dset_handle = H5Dcreate (m_handle, name.c_str(), memtype, dspace_handle, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            status = H5Dwrite(dset_handle, memtype, dspace_handle, dspace_handle, H5P_DEFAULT, buffer.data());
            REQUIRE_FALSE_ALL(status, "HDF5 Error on string array save");
            H5Dclose(dset_handle);
            H5Sclose(dspace_handle);
        }

        /**
         * wrapper for save in the case that only a single string is to be stored
         */
        void save(const str_t& name, const str_t& v, uint_t irank=0ul) const {
            strv_t vs = {v};
            save(name, vs, irank);
        }
    };

}


#endif //M7_HDF5_NODE_H
