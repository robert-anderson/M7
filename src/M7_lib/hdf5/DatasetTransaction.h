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
        hid_t m_plist = 0;
        hid_t m_filespace = 0;
        hid_t m_file_hyperslab = 0;
        hid_t m_mem_hyperslab = 0;
        hid_t m_dataset = 0;
    public:
        DatasetTransaction(dataset::DistListFormat format):
            m_format(std::move(format)), m_counts(m_format.m_local.ndim(), 0), m_offsets(m_format.m_local.ndim(), 0){}

        ~DatasetTransaction() {
            // free all HDF5 handles if they have been opened (set non-zero) by the subclasses
            if (m_plist) H5Pclose(m_plist);
            else logging::warn("property list was not opened in transaction");
            if(m_filespace) H5Sclose(m_filespace);
            else logging::warn("filespace was not opened in transaction");
            if (m_file_hyperslab) H5Sclose(m_file_hyperslab);
            else logging::warn("file hyperslab was not opened in transaction");
            if (m_mem_hyperslab) H5Sclose(m_mem_hyperslab);
            else logging::warn("application memory hyperslab was not opened in transaction");
            if (m_dataset) H5Dclose(m_dataset);
            else logging::warn("dataset was not opened in transaction");
        }

        /*
         * since we have a non-default dtor, follow the rule of three...
         */
        DatasetTransaction(const DatasetTransaction& other): DatasetTransaction(other.m_format){}
        DatasetTransaction& operator=(const DatasetTransaction&) {
            return *this;
        }

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

        /**
         * to keep down on the number of optional arguments to the save_* methods
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
        static void save_dist_list(
            const NodeWriter& nw,
            const str_t& name,
            const void* src,
            Type type,
            bool is_complex,
            uintv_t item_shape,
            strv_t item_dim_names,
            uint_t nitem,
            const Options& opts);

        /**
         * type-templated wrapper of the above, raw-data function
         */
        template<typename T>
        static void save_dist_list(
            const NodeWriter& nw,
            const str_t& name,
            const T* src,
            uintv_t item_shape,
            strv_t item_dim_names,
            uint_t nitem,
            const Options& opts)
        {
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem, opts);
        }

        /**
         * infer local nitem based on the length of the source vector and whether this MPI rank is involved in the save
         */
        template<typename T>
        static void save_dist_list(
            const NodeWriter& nw,
            const str_t& name,
            const v_t<T>& src,
            uintv_t item_shape,
            strv_t item_dim_names,
            bool this_rank,
            const Options& opts)
        {
            const auto nitem = this_rank ? src.size() / nd::nelement(item_shape) : 0ul;
            REQUIRE_FALSE(src.size() % nd::nelement(item_shape), "item shape inconsistent with size of data buffer");
            save_dist_list(nw, name, src.data(), item_shape, item_dim_names, nitem, opts);
        }

        /**
         * src is not a distributed list, so reinterpret the leading dimension as a distributed dimension and delegate
         * the save_dist_list method
         */
        template<typename T>
        static void save_array(
            const NodeWriter& nw,
            const str_t& name,
            const T* src,
            uintv_t shape,
            strv_t dim_names,
            uint_t irank=0ul,
            const Options& opts={})
        {
            auto nitem = mpi::i_am(irank) ? shape.front() : 0ul;
            uintv_t item_shape(shape.cbegin()+1, shape.cend());
            auto leading_dim_name = dim_names.empty() ? "" : dim_names.front();
            strv_t item_dim_names(dim_names.cbegin()+1, dim_names.cend());
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem, opts);
        }

        template<typename T>
        static void save_array(
            const NodeWriter& nw,
            const str_t& name,
            const v_t<T>& src,
            uintv_t shape,
            strv_t dim_names,
            uint_t irank=0ul,
            const Options& opts={})
        {
            if (mpi::i_am(irank)) {
                REQUIRE_EQ(nd::nelement(shape), src.size(), "number of elements in array is incompatible with given shape");
            }
            save_array(nw, name, src.data(), shape, dim_names, irank, opts);
        }

        template<typename T>
        static void save_vector(
            const NodeWriter& nw,
            const str_t& name,
            const v_t<T>& src,
            uint_t irank=0ul,
            const Options& opts={})
        {
            save_array(nw, name, src, {src.size()}, {"elements"}, irank, opts);
        }

        template<typename T>
        static void save_scalar(const NodeWriter& nw, const str_t& name, const T& src, uint_t irank=0ul, const Options& opts={}){
            save_array(nw, name, &src, {1ul}, {"scalar"}, irank, opts);
        }
    };

    class DatasetLoader : public DatasetTransaction {

        static bool exists(hid_t parent, const str_t &name);

        static void require_exists(hid_t parent, const str_t &name);

        static uint_t read_dataset_ndim(hid_t parent, const str_t &name);

        static uintv_t read_shape(hid_t parent, const str_t &name);

        static hdf5::dataset::FullDistListFormat read_full_format(hid_t parent, const str_t &name, bool this_rank);

        static hdf5::dataset::PartDistListFormat read_part_format(hid_t parent, const str_t &name, bool this_rank);

        static bool valid_part_flag(bool part);

        /**
         * @param parent
         *  Node (group or file) open for reading
         * @param name
         *  name of the dataset whose sizes are queried
         * @param part
         *  true if this is a partial read, where each rank with this_rank=true has a portion of the data to read
         *  false if every rank with this_rank=true reads the entire dataset
         * @param this_rank
         *  true if this MPI rank participates in the load operation
         * @return
         *  pair of integers specifying the size of the locally-read data: nitem and item_size
         */
        static dataset::DistListFormat read_format(hid_t parent, const str_t &name, bool part, bool this_rank);

    public:

        DatasetLoader(const NodeReader &nr, const str_t &name, bool part, bool this_rank);

        void load_attrs(std::list<Attr>& attrs) const;

        bool read(void* dst, uint_t nitem);

        ~DatasetLoader();

        /**
         * to keep down on the number of optional arguments to the load_* methods
         */
        struct Options {
            uint_t m_max_nitem_per_op;
            bool m_part;
            bool m_this_rank;
            Options(bool part, bool this_rank):
                m_max_nitem_per_op(c_default_max_nitem_per_op), m_part(part), m_this_rank(this_rank){}
        };

        /**
         * load data from hdf5 dataset into distributed, contiguous raw buffers
         * @param nr
         *  node reader object of the dataset's parent
         * @param name
         *  name in the HDF5 node at which the dataset is stored
         * @param dst
         *  pointer to the beginning of the buffer into which the entire dataset (or this rank's share thereof) is to
         *  be read
         * @param size
         *  size of buffer in bytes (must be at least as large as that required by the dataset)
         * @param type
         *  native type of data stored at dst
         * @param is_complex
         *  true if the data stored at dst is of std::complex<T> type
         * @param opts
         *  passed as mutable reference so that
         * @param attrs
         *  all attributes linked to the dataset
         */
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            void* dst,
            uint_t size,
            hdf5::Type type,
            bool is_complex,
            const Options& opts,
            std::list<Attr>& attrs);

        /**
         * type-templated wrapper of the above, raw-data function
         */
        template<typename T>
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            T* dst,
            uint_t size,
            const Options& opts,
            std::list<Attr>& attrs)
        {
            load_dist_list(nr, name, dst, size*sizeof(T), Type::make<T>(), dtype::is_complex<T>(), opts, attrs);
        }

        /**
         * as above, but discarding attrs
         */
        template<typename T>
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            T* dst,
            uint_t size,
            const Options& opts)
        {
            std::list<Attr> attrs;
            load_dist_list(nr, name, dst, size, opts, attrs);
        }

        /**
         * infer size of vector from dataset
         */
        template<typename T>
        static void load_dist_list(
                const NodeReader& nr, const str_t& name, v_t<T>& dst, const Options& opts, std::list<Attr>& attrs){
            const auto size = read_format(nr.m_handle, name, opts.m_part, opts.m_this_rank).m_local.m_size;
            dst.resize(size / sizeof(T));
            load_dist_list(nr, name, dst, dst.size(), opts, attrs);
        }

        /**
         * as above, but discarding attrs
         */
        template<typename T>
        static void load_dist_list(const NodeReader& nr, const str_t& name, v_t<T>& dst, const Options& opts){
            std::list<Attr> attrs;
            load_dist_list(nr, name, dst, opts, attrs);
        }
    };
}


#endif //M7_DATASETTRANSACTION_H
