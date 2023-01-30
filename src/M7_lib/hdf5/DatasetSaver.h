//
// Created by rja on 30/01/23.
//

#ifndef M7_DATASETSAVER_H
#define M7_DATASETSAVER_H

#include <utility>

#include "NodeWriter.h"
#include "M7_lib/util/Vector.h"

namespace hdf5 {
    struct DatasetSaver {
        const hdf5::dataset::PartDistListFormat m_format;

    private:
        /**
         * counter for the number of items already transacted
         */
        uint_t m_nitem_saved = 0ul;
        /**
         * working vectors for determining hyperslabs
         */
        v_t<hsize_t> m_counts, m_offsets;

        /**
         * HDF5 asset handles
         */
        hid_t m_plist;
        hid_t m_filespace;
        hid_t m_file_hyperslab;
        hid_t m_mem_hyperslab;
        hid_t m_dataset;
    public:

        DatasetSaver(const NodeWriter& nw, const str_t& name, hdf5::dataset::PartDistListFormat format):
        m_format(std::move(format)){
            REQUIRE_FALSE_ALL(name.empty(), "HDF5 dataset must be given a name")
            m_filespace = H5Screate_simple(format.m_h5_shape.size(), format.m_h5_shape.data(), nullptr);

            // specify format of the dataset
            m_dataset = H5Dcreate(nw.m_handle, name.c_str(), format.m_local.m_item.m_type, m_filespace,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            REQUIRE_GT_ALL(m_dataset, 0, "dataset creation failed");
            {
                // set dimension labels
                uint_t i = 0ul;
                for (auto& dim_name: format.m_local.m_dim_names) H5DSset_label(m_dataset, i++, dim_name.c_str());
            }
            H5Sclose(m_filespace);
            m_filespace = H5Dget_space(m_dataset);

            // datasets are always transacted in collective I/O mode
            m_plist = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(m_plist, H5FD_MPIO_COLLECTIVE);

            const hsize_t ndim = format.m_h5_shape.size();
            // initialize an array to specify the counts of data transacted. this only changes in the first element
            m_counts = vector::prepended(format.m_local.m_item.m_h5_shape, 0ul);
            // hyperslabs for the offsets and counts associated with the file and memory layout respectively
            m_file_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);
            m_mem_hyperslab = H5Screate_simple(ndim, format.m_h5_shape.data(), nullptr);

            // first element of offset is incremented at each iteration
            m_offsets = v_t<hsize_t>(ndim);

            if (!m_format.m_nitem) logging::warn("Saving empty dataset \"{}\"", name);
        }

        /**
         * save distributed data formatted as contiguous raw buffers in one write operation
         * @param nw
         *  node writer object of the dataset's parent
         * @param name
         *  name to be given to the dataset in the node writer
         * @param src
         *  local buffer to be written
         * @param type
         *  HDF5 native type of the numbers represented in src
         * @param is_complex
         *  if true, the src is taken to be casted from a std::complex<T>, where T determines the above type
         * @param item_shape
         *  shape of each item
         * @param item_dim_names
         *  labels associated with each item dimension (does not include the leading dimension or real/imag)
         * @param nitem
         *  number of items stored on this MPI rank
         * @param attrs
         *  attributes to save inside the dataset HDF5 node
         * @param leading_dim_name
         *  name to give to the distributed dimension
         */
        static void save_dist_list(const NodeWriter& nw, const str_t& name, const void* src, Type type, bool is_complex,
                         uintv_t item_shape, strv_t item_dim_names, uint_t nitem,
                         std::list<Attr> attrs = {}, str_t leading_dim_name="item") {
            const dataset::ItemFormat item_format(type, std::move(item_shape), std::move(item_dim_names), is_complex);
            const dataset::PartDistListFormat format(item_format, nitem, std::move(leading_dim_name));
            DatasetSaver ds(nw, name, format);
            ds.write(src, nitem);
            ds.save_attrs(attrs);
        }

        /**
         * type-templated wrapper of the above, raw-data function
         */
        template<typename T>
        static void save_dist_list(const NodeWriter& nw, const str_t& name, const T* src, uintv_t item_shape,
                         strv_t item_dim_names, uint_t nitem, std::list<Attr> attrs = {}, str_t leading_dim_name="item"){
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(),
                    item_shape, item_dim_names, nitem, attrs, leading_dim_name);
        }

        /**
         * infer local nitem based on the length of the source vector and whether this MPI rank is involved in the save
         */
        template<typename T>
        static void save_dist_list(const NodeWriter& nw, const str_t& name, const v_t<T>& src, uintv_t item_shape,
                                   strv_t item_dim_names, bool this_rank, std::list<Attr> attrs = {}, str_t leading_dim_name="item"){
            const auto nitem = this_rank ? src.size() / nd::nelement(item_shape) : 0ul;
            REQUIRE_FALSE(src.size() % nd::nelement(item_shape), "item shape inconsistent with size of data buffer");
            save_dist_list(nw, name, src.data(), item_shape, item_dim_names, nitem, attrs, leading_dim_name);
        }

        /**
         * src is not a distributed list, so reinterpret the leading dimension as a distributed dimension and delegate
         * the save_dist_list method
         */
        template<typename T>
        static void save_array(const NodeWriter& nw, const str_t& name, const T* src, uintv_t shape,
                               strv_t dim_names, uint_t irank=0ul, std::list<Attr> attrs = {}){
            auto nitem = mpi::i_am(irank) ? shape.front() : 0ul;
            uintv_t item_shape(shape.cbegin()+1, shape.cend());
            auto leading_dim_name = dim_names.empty() ? "" : dim_names.front();
            strv_t item_dim_names(dim_names.cbegin()+1, dim_names.cend());
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem, attrs, leading_dim_name);
        }


        /**
         * write raw data collectively
         * @param src
         *  pointer to contiguous raw data buffer
         * @param nitem
         *  number of items stored in raw data
         * @return
         *  true if all ranks are done writing
         */
        bool write(const void* src, uint_t nitem) {
            /*
             * only call H5Dwrite if there is any data across all MPI ranks
             */
            if (!m_format.m_nitem) return true;
            // the number of items transacted may not be larger than the number of items remaining
            REQUIRE_LE(m_nitem_saved + nitem, m_format.m_local.m_nitem, "too many items");
            m_counts[0] = nitem;
            m_offsets[0] = 0;
            H5Sselect_hyperslab(m_mem_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
            m_offsets[0] = m_format.m_nitem_displ + m_nitem_saved;
            REQUIRE_EQ(bool(src), bool(m_counts[0]), "count zero with non-null data or count non-zero with null data");
            H5Sselect_hyperslab(m_file_hyperslab, H5S_SELECT_SET, m_offsets.data(), nullptr, m_counts.data(), nullptr);
            auto status = H5Dwrite(m_dataset, m_format.m_local.m_item.m_type, m_mem_hyperslab, m_file_hyperslab, m_plist, src);
            REQUIRE_FALSE(status, "HDF5 Error on multidimensional save");
            bool all_done = !src;
            all_done = mpi::all_land(all_done);
            // deplete number of remaining items by the number just transacted
            m_nitem_saved += m_counts[0];
            return all_done;
        }

        void save_attrs(const std::list<Attr>& attrs) {
            for (const auto& attr: attrs) attr.save(m_dataset);
        }

        ~DatasetSaver() {
            // free all HDF5 handles
            H5Pclose(m_plist);
            H5Sclose(m_filespace);
            H5Sclose(m_file_hyperslab);
            H5Sclose(m_mem_hyperslab);
            H5Dclose(m_dataset);
        }
    };
}


#endif //M7_DATASETSAVER_H
