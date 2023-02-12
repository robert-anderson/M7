//
// Created by rja on 04/02/23.
//

#ifndef M7_DATASETTRANSACTION_H
#define M7_DATASETTRANSACTION_H

#include "DatasetFormat.h"
#include "Node.h"

namespace hdf5 {
    static constexpr uint_t c_default_max_nitem_per_op = 16000000;
    static constexpr char c_default_leading_dim_name[] = "item";
    /**
     * common members between dataset saver and loader
     */
    struct DatasetTransaction {
        const hdf5::dataset::DistListFormat m_format;
    protected:
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
        explicit DatasetTransaction(dataset::DistListFormat format);

        ~DatasetTransaction();

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
         * @param max_nitem_per_op
         *  maximum number of items to transact in a single write call
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
            std::list<Attr> attrs = {},
            uint_t max_nitem_per_op = c_default_max_nitem_per_op,
            str_t leading_dim_name = c_default_leading_dim_name);

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
            std::list<Attr> attrs = {},
            uint_t max_nitem_per_op = c_default_max_nitem_per_op,
            str_t leading_dim_name = c_default_leading_dim_name)
        {
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem,
                           attrs, max_nitem_per_op, leading_dim_name);
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
            std::list<Attr> attrs = {},
            uint_t max_nitem_per_op = c_default_max_nitem_per_op,
            str_t leading_dim_name = c_default_leading_dim_name)
        {
            const auto nitem = this_rank ? src.size() / nd::nelement(item_shape) : 0ul;
            REQUIRE_FALSE(src.size() % nd::nelement(item_shape), "item shape inconsistent with size of data buffer");
            save_dist_list(nw, name, src.data(), item_shape, item_dim_names, nitem,
                           attrs, max_nitem_per_op, leading_dim_name);
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
            std::list<Attr> attrs = {},
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            auto nitem = mpi::i_am(irank) ? shape.front() : 0ul;
            uintv_t item_shape(shape.cbegin()+1, shape.cend());
            auto leading_dim_name = dim_names.empty() ? "" : dim_names.front();
            strv_t item_dim_names(dim_names.cbegin()+1, dim_names.cend());
            save_dist_list(nw, name, src, Type::make<T>(), dtype::is_complex<T>(), item_shape, item_dim_names, nitem,
                           attrs, max_nitem_per_op);
        }

        template<typename T>
        static void save_array(
            const NodeWriter& nw,
            const str_t& name,
            const v_t<T>& src,
            uintv_t shape,
            strv_t dim_names,
            uint_t irank=0ul,
            std::list<Attr> attrs = {},
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            if (mpi::i_am(irank)) {
                REQUIRE_EQ(nd::nelement(shape), src.size(), "number of elements in array is incompatible with given shape");
            }
            save_array(nw, name, src.data(), shape, dim_names, irank, attrs, max_nitem_per_op);
        }

        template<typename T>
        static void save_vector(
            const NodeWriter& nw,
            const str_t& name,
            const v_t<T>& src,
            uint_t irank=0ul,
            std::list<Attr> attrs = {},
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            save_array(nw, name, src, {src.size()}, {"elements"}, irank, attrs, max_nitem_per_op);
        }

        template<typename T>
        static void save_scalar(const NodeWriter& nw, const str_t& name, const T& src, uint_t irank=0ul){
            save_array(nw, name, &src, {1ul}, {"scalar"}, irank);
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

    public:

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

        DatasetLoader(const NodeReader &nr, const str_t &name, bool part, bool this_rank);

        void load_attrs(std::list<Attr>& attrs) const;

        bool read(void* dst, uint_t nitem);

        ~DatasetLoader();

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
         * @param part
         *  true if the loading is to be split up equally among MPI ranks, false if all ranks read all data
         * @param this_rank
         *  true if this MPI rank participates in the load operation (reads partial or full dataset from file)
         * @param attrs
         *  all attributes linked to the dataset
         * @param max_nitem_per_op
         *  maximum number of items to transfer in a single read operation
         */
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            void* dst,
            uint_t size,
            hdf5::Type type,
            bool is_complex,
            bool part,
            bool this_rank,
            std::list<Attr>& attrs,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op);

        /**
         * type-templated wrapper of the above, raw-data function
         */
        template<typename T>
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            T* dst,
            uint_t size,
            bool part,
            bool this_rank,
            std::list<Attr>& attrs,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            load_dist_list(nr, name, dst, size*sizeof(T), Type::make<T>(), dtype::is_complex<T>(),
                part, this_rank, attrs, max_nitem_per_op);
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
            bool part,
            bool this_rank,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            std::list<Attr> attrs;
            load_dist_list(nr, name, dst, size, part, this_rank, attrs, max_nitem_per_op);
        }

        /**
         * infer size of vector from dataset
         */
        template<typename T>
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            v_t<T>& dst,
            bool part,
            bool this_rank,
            std::list<Attr>& attrs,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            const auto size = read_format(nr.m_handle, name, part, this_rank).m_local.m_size;
            dst.resize(size / sizeof(T));
            load_dist_list(nr, name, dst.data(), dst.size(), part, this_rank, attrs, max_nitem_per_op);
        }

        /**
         * as above, but discarding attrs
         */
        template<typename T>
        static void load_dist_list(
            const NodeReader& nr,
            const str_t& name,
            v_t<T>& dst,
            bool part,
            bool this_rank,
            uint_t max_nitem_per_op = c_default_max_nitem_per_op)
        {
            std::list<Attr> attrs;
            load_dist_list(nr, name, dst, part, this_rank, attrs, max_nitem_per_op);
        }

        /**
         * convenience function for reading into and returning a std::vector
         */
        template<typename T>
        static v_t<T> load_vector(const NodeReader& nr, const str_t& name, bool this_rank=true){
            v_t<T> tmp;
            load_dist_list(nr, name, tmp, false, this_rank);
            return tmp;
        }
    };
}


#endif //M7_DATASETTRANSACTION_H
