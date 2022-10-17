//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_NDDISTLIST_H
#define M7_HDF5_NDDISTLIST_H


#include "File.h"


namespace hdf5 {
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
        const v_t<hsize_t> m_item_dims;
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
        const v_t<hsize_t> m_list_dims_local;
        /**
         * global shape of the list
         */
        const v_t<hsize_t> m_list_dims_global;
        /**
         * sum of m_nitem_local for all MPI ranks with index less than this rank
         */
        const hsize_t m_item_offset;

        /**
         * extent of the currently selected hyperslab in the HDF5 dataset
         */
        v_t<hsize_t> m_hyperslab_counts;
        /**
         * multidimensional offset for the currently selected hyperslab
         */
        v_t<hsize_t> m_hyperslab_offsets;
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
         * dtype as an integer
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
        v_t<hsize_t> get_list_dims_local();

        /**
         * stick the global number of items onto the front of the item shape
         * @return
         *  list dims aka the overall shape of the global dataset
         */
        v_t<hsize_t> get_list_dims_global();

        /**
         * @return
         *  flat offset to the first item handled by this rank
         */
        hsize_t get_item_offset();

        NdDistListBase(hid_t parent_handle, str_t name, const uintv_t &item_dims, const uint_t &nitem,
                       bool writemode, hid_t h5type);

        /**
         * sets the internal state to select the iitem-th item on this rank
         * @param iitem
         *  item index for selection
         */
        void select_hyperslab(const uint_t &iitem);

        ~NdDistListBase();
    };

    /**
     * The writing case for a distributed N-dimensional list
     */
    struct NdDistListWriter : public NdDistListBase {
    private:
        NdDistListWriter(hid_t parent_handle, str_t name, const uintv_t &item_dims,
                         const uint_t &nitem, hid_t h5type, const strv_t &dim_labels = {});

    public:
        NdDistListWriter(NodeWriter &parent, str_t name, const uintv_t &item_dims,
                         const uint_t &nitem, hid_t h5type, const strv_t &dim_labels = {}) :
                NdDistListWriter(parent.m_handle, name, item_dims, nitem, h5type, dim_labels) {}

        /**
         * write the item stored at data to the iitem position within this rank
         * @param iitem
         *  item index determining the position on the HDF5 file buffer into which the data is copied
         * @param data
         *  pointer to data of any type. the counts in the hyperslab selection determine the number of elements n of the
         *  native type T corresponding to m_h5type are to be written. Thus the sizeof(T)*n bytes after data are copied
         */
        void write_h5item_bytes(const uint_t &iitem, const void *data);

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
        static uint_t extract_list_ndim(hid_t parent_handle, str_t name) {
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
        static uintv_t extract_list_dims(hid_t parent_handle, str_t name) {
            auto rank = extract_list_ndim(parent_handle, name);
            auto dataset = H5Dopen1(parent_handle, name.c_str());
            auto dataspace = H5Dget_space(dataset);
            v_t<hsize_t> dims(rank, 0ul);
            H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            uintv_t out;
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
        static uintv_t extract_item_dims(hid_t parent_handle, str_t name) {
            auto list_dims = extract_list_dims(parent_handle, name);
            uintv_t out;
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
        static uint_t extract_nitem(hid_t parent_handle, str_t name) {
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
        static uint_t local_nitem(hid_t parent_handle, str_t name) {
            auto nitem_tot = extract_nitem(parent_handle, name);
            auto nitem_share = nitem_tot / mpi::nrank();
            // give the remainder to the root rank
            if (mpi::i_am_root()) return nitem_share + nitem_tot % mpi::nrank();
            else return nitem_share;
        }

    private:
        NdDistListReader(hid_t parent_handle, str_t name, hid_t h5_type) :
                NdDistListBase(parent_handle, name, extract_item_dims(parent_handle, name),
                               local_nitem(parent_handle, name), false, h5_type) {}

    public:
        NdDistListReader(NodeReader &parent, str_t name, hid_t h5_type) :
                NdDistListReader(parent.m_handle, name, h5_type) {}

        /**
         * read the item stored at the iitem position into the data pointer
         * @param iitem
         *  item index determining the position on the HDF5 file buffer from which the data is copied
         * @param data
         *  pointer to data of any type. the counts in the hyperslab selection determine the number of elements n of the
         *  native type T corresponding to m_h5type are to be read.
         */
        void read_h5item_bytes(const uint_t &iitem, void *data) {
            select_hyperslab(iitem);
            logging::debug_("reading data...");
            if (data) {
                auto status = H5Dread(m_dataset_handle, m_h5type, m_memspace_handle,
                                      m_filespace_handle, m_coll_plist, data);
                DEBUG_ONLY(status);
                DEBUG_ASSERT_TRUE(!status, "HDF5 read failed");
            } else {
                auto status = H5Dread(m_dataset_handle, m_h5type, m_none_memspace_handle,
                                      m_filespace_handle, m_coll_plist, data);
                DEBUG_ONLY(status);
                DEBUG_ASSERT_TRUE(!status, "HDF5 read failed");
            }
            logging::debug_("data read");
        }
    };
}

#endif //M7_HDF5_NDDISTLIST_H
