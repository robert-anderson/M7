//
// Created by rja on 04/02/23.
//

#ifndef M7_DATASETTRANSACTION_H
#define M7_DATASETTRANSACTION_H

#include "DatasetFormat.h"
#include "NodeWriter.h"
#include "NodeReader.h"

namespace hdf5 {
    static constexpr uint_t c_default_max_nitem_per_op = 16000000;
    /**
     * common members between dataset saver and loader
     */
    struct DatasetTransaction {
    protected:
        const hdf5::dataset::DistListFormat m_format;
        /**
         * counter for the number of items already transacted
         */
        uint_t m_nitem_done = 0ul;
        /**
         * working vectors for determining hyperslabs
         */
        v_t<hsize_t> m_counts, m_offsets;

        /**
         * HDF5 object handles
         */
        hid_t m_plist;
        hid_t m_filespace;
        hid_t m_file_hyperslab;
        hid_t m_mem_hyperslab;
        hid_t m_dataset;
    public:
        DatasetTransaction(dataset::DistListFormat format):
            m_format(std::move(format)), m_counts(m_format.m_local.ndim(), 0), m_offsets(m_format.m_local.ndim(), 0){}


        uint_t nitem_done() const {
            return m_nitem_done;
        }

        uint_t nitem_remaining() const {
            return m_format.m_local.m_nitem - m_nitem_done;
        }

        uint_t nitem_next(uint_t max_nitem_per_op) const {
            return std::min(nitem_remaining(), max_nitem_per_op);
        }
    };

    struct DatasetSaver : public DatasetTransaction {
        DatasetSaver(const NodeWriter& nw, const str_t& name, hdf5::dataset::PartDistListFormat format);

        /**
         * to keep cut down on the number of optional arguments to the save_* methods
         */
        struct Options {
            uint_t m_max_nitem_per_op;
            std::list<Attr> m_attrs;
            str_t m_leading_dim_name;
            Options(): m_max_nitem_per_op(c_default_max_nitem_per_op), m_attrs(), m_leading_dim_name("item"){}
        };

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
                                   uintv_t item_shape, strv_t item_dim_names, uint_t nitem, const Options& opts);

        /**
         * type-templated wrapper of the above, raw-data function
         */
        template<typename T>
        static void save_dist_list(const NodeWriter& nw, const str_t& name, const T* src, uintv_t item_shape,
                                   strv_t item_dim_names, uint_t nitem, const Options& opts){
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem, opts);
        }

        /**
         * infer local nitem based on the length of the source vector and whether this MPI rank is involved in the save
         */
        template<typename T>
        static void save_dist_list(const NodeWriter& nw, const str_t& name, const v_t<T>& src, uintv_t item_shape,
                                   strv_t item_dim_names, bool this_rank, const Options& opts){
            const auto nitem = this_rank ? src.size() / nd::nelement(item_shape) : 0ul;
            REQUIRE_FALSE(src.size() % nd::nelement(item_shape), "item shape inconsistent with size of data buffer");
            save_dist_list(nw, name, src.data(), item_shape, item_dim_names, nitem, opts);
        }

        /**
         * src is not a distributed list, so reinterpret the leading dimension as a distributed dimension and delegate
         * the save_dist_list method
         */
        template<typename T>
        static void save_array(const NodeWriter& nw, const str_t& name, const T* src, uintv_t shape, strv_t dim_names,
                               uint_t irank=0ul, const Options& opts={}){
            auto nitem = mpi::i_am(irank) ? shape.front() : 0ul;
            uintv_t item_shape(shape.cbegin()+1, shape.cend());
            auto leading_dim_name = dim_names.empty() ? "" : dim_names.front();
            strv_t item_dim_names(dim_names.cbegin()+1, dim_names.cend());
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem, opts);
        }

        template<typename T>
        static void save_array(const NodeWriter& nw, const str_t& name, const v_t<T>& src, uintv_t shape,
                               strv_t dim_names, uint_t irank=0ul, const Options& opts={}){
            if (mpi::i_am(irank)) {
                REQUIRE_EQ(nd::nelement(shape), src.size(), "number of elements in array is incompatible with given shape");
            }
            save_array(nw, name, src.data(), shape, dim_names, irank, opts);
        }

        template<typename T>
        static void save_vector(const NodeWriter& nw, const str_t& name, const v_t<T>& src,
                                uint_t irank=0ul, const Options& opts={}){
            save_array(nw, name, src, {src.size()}, {"elements"}, irank, opts);
        }

        template<typename T>
        static void save_scalar(const NodeWriter& nw, const str_t& name, const T& src,
                                uint_t irank=0ul, const Options& opts={}){
            save_array(nw, name, &src, {1ul}, {"scalar"}, irank, opts);
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
        bool write(const void* src, uint_t nitem);

        void save_attrs(const std::list<Attr>& attrs);

        ~DatasetSaver();
    };
};


#endif //M7_DATASETTRANSACTION_H
