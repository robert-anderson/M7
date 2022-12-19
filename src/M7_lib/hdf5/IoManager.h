//
// Created by rja on 06/12/22.
//

#ifndef M7_HDF5_BUFFERER_H
#define M7_HDF5_BUFFERER_H

#include <utility>

#include "Type.h"
#include "M7_lib/util/Vector.h"
#include "M7_lib/util/Pointer.h"

/**
 *
 */
namespace hdf5 {

    /**
     * we require a set of classes to abstract the process of dataset I/O
     *
     * write data from a single designated MPI rank callable:
     *  - collectively
     *  - independently (on the definitive rank only)
     * write data from all ranks into the same dataset collectively
     * read all data to a single designated rank
     *  - collectively
     *  - independently (on the designated rank only)
     * read data to all ranks
     */

//    enum IoPolicy {
//        WriteFromAll,
//        WriteFromOneColl,
//        WriteFromOneIndep,
//        ReadToAll,
//        ReadToOneColl,
//        ReadToOneIndep,
//    };

    namespace dataset {

        /**
         * nomenclature:
         *  type = an HDF5 native type (cannot be complex)
         *  elem = single number as stored in the C++ program (real or complex)
         *  h5_elem = single number as stored in the HDF5 archive (cannot be complex)
         *  item = a shaped array of elems
         *  h5_item = a shaped array of h5_elems
         */
        struct Format {
            /**
             * native type of the HDF5 dataset involved in the I/O
             */
            const Type m_h5_type;
            /**
             * layout of an item in C++ terms
             */
            const uintv_t m_item_shape;
            /**
             * total number of items (locally held)
             */
            const uint_t m_nitem;
            /**
             * layout of an item in HDF5 terms in hsize_t type
             */
            const v_t<hsize_t> m_h5_item_shape;
            /**
             * size in bytes of each item
             */
            const uint_t m_item_size;
            /**
             * layout of all locally stored items in HDF5 terms in hsize_t type
             */
            const v_t<hsize_t> m_h5_shape;
            /**
             * names of each dimension in the HDF5 metadata
             */
            const strv_t m_h5_dim_names;

            Format(Type h5_type, uintv_t item_shape, uint_t nitem, strv_t dim_names, bool add_complex_dim) :
                m_h5_type(h5_type), m_item_shape(std::move(item_shape)), m_nitem(nitem),
                m_h5_item_shape(
                        convert::vector<hsize_t>(add_complex_dim ? vector::appended(m_item_shape, 2ul) : m_item_shape)),
                m_item_size(nd::nelement(m_h5_item_shape) * m_h5_type.m_size),
                m_h5_shape(convert::vector<hsize_t>(vector::prepended(m_h5_item_shape, m_nitem))),
                m_h5_dim_names(add_complex_dim && !dim_names.empty() ? vector::appended(dim_names, "real/imag") : dim_names) {
                REQUIRE_TRUE(m_h5_dim_names.empty() || m_h5_dim_names.size() == m_h5_item_shape.size(),
                             "incorrect number of dimension names");
                REQUIRE_TRUE((m_h5_item_shape.back()==2ul) || !add_complex_dim,
                             "last item in item shape must be 2 if data is add_complex_dim-valued");
            }
        };

        struct DistFormat : Format{
            const uint_t m_nitem_sum;
            const uintv_t m_nitem_all;
            const uintv_t m_nitem_offsets;
            const v_t<hsize_t> m_h5_shape_sum;
            DistFormat(Type h5_type, uintv_t item_shape, uint_t nitem, strv_t dim_names, bool add_complex_dim):
                Format(h5_type, std::move(item_shape), nitem, std::move(dim_names), add_complex_dim),
                m_nitem_sum(mpi::all_sum(m_nitem)), m_nitem_all(mpi::all_gathered(m_nitem)),
                m_nitem_offsets(mpi::counts_to_displs_consec(m_nitem_all)),
                m_h5_shape_sum(vector::prepended(m_h5_item_shape, m_nitem_sum)){}
        };

        typedef std::function<const void*(uint_t i, const dataset::DistFormat& format, uint_t max_nitem_per_op)> save_fn;
        typedef std::function<void*(uint_t i, const dataset::DistFormat& format, uint_t max_nitem_per_op)> load_fn;

        struct SaveManager {
            const DistFormat m_format;
            const uint_t m_max_nitem_per_write;
            SaveManager(DistFormat format, uint_t max_nitem_per_write):
                m_format(std::move(format)), m_max_nitem_per_write(max_nitem_per_write){}
            /**
             * prepare the buffer for the next write operation
             * @return
             *  pointer to the bytes to be written, null if there is no more data on this MPI rank
             */
            virtual const void* get(uint_t iblock) const = 0;
        };
        template<typename T>
        struct ContiguousSaveManager : SaveManager {
            const char* const m_begin;
            const char* const m_end;
            ContiguousSaveManager(const T* v, uint_t size, uint_t max_nitem_per_op): SaveManager(
                {hdf5::Type::make<T>(), {1ul}, size, {"element"}, dtype::is_complex<T>()}, max_nitem_per_op),
                m_begin(reinterpret_cast<const char*>(v)), m_end(m_begin + m_format.m_item_size * m_format.m_nitem){}

            ContiguousSaveManager(const T* v, uint_t size): ContiguousSaveManager(v, size, size){}

            ContiguousSaveManager(const v_t<T>* v): ContiguousSaveManager(v->data(), v->size()){}

            const void* get(uint_t i) const override {
                const auto ptr = m_begin + (i * m_max_nitem_per_write * m_format.m_item_size);
                return ::ptr::in_range(ptr, m_begin, m_end) ? ptr : nullptr;
            }
        };

        template<typename T>
        ContiguousSaveManager<T> make_saver(const T* v, uint_t size, uint_t max_nitem_per_op) {
            return {v, size, max_nitem_per_op};
        }

        template<typename T>
        ContiguousSaveManager<T> make_saver(const T* v, uint_t size) {
            return make_saver(v, size, size);
        }

        template<typename T>
        ContiguousSaveManager<T> make_saver(const v_t<T>* v) {
            return make_saver(v, v->size());
        }

//        struct LoadManager : IoManager {
//            LoadManager(Format format, uint_t max_nitem_per_op): IoManager(std::move(format), max_nitem_per_op){}
//            /**
//             * prepare the buffer for the next write operation
//             * @return
//             *  pointer to the bytes to be written, null if there is no more data on this MPI rank
//             */
//            virtual void* get_next() = 0;
//        };

//        template<typename T>
//        struct VectorLoadManager : LoadManager {
//            const void* const m_begin;
//            const void* const m_end;
//            VectorLoadManager(v_t<T>* v, uint_t max_nitem_per_op): SaveManager(
//                {hdf5::Type::make<T>(), {1ul}, v->size(), {"element"}, dtype::is_complex<T>()}, max_nitem_per_op),
//                m_begin(reinterpret_cast<void const*>(v->data())),
//                m_end(m_begin+m_format.m_item_size*m_format.m_nitem){}
//
//            VectorLoadManager(const v_t<T>* v): VectorSaveManager(v, v->size()){}
//
//            const void * get_next() override {
//                auto ptr = m_begin;
//                ptr += m_nop*m_max_nitem_per_op*m_format.m_item_size;
//                return ::ptr::in_range(ptr, m_begin, m_end) ? ptr : nullptr;
//            }
//        };

//        struct DistFormat : Format {
//            v_t
//            DistFormat(Type h5_type, uintv_t item_shape, uint_t nitem, strv_t dim_names, bool complex) :
//                Format(h5_type, std::move(item_shape), nitem, std::move(dim_names), complex){}
//        };

    }

#if 0

    struct WriteManager: IoManager {
        WriteManager(Type type, uintv_t item_shape, uint_t nitem, uint_t nitem_per_transfer, strv_t dim_names, bool complex):
                IoManager(type, item_shape, nitem, nitem_per_transfer, dim_names, complex){}

        virtual void* transfer(uint_t& nitem) = 0;
    };

    template<typename T>
    struct VectorWriteManager : WriteManager {
        const v_t<T>* m_vec;
        VectorWriteManager(const v_t<T>* vec, uint_t nitem_per_transfer):
                WriteManager(Type(vec->data()), {}, vec->size(), nitem_per_transfer, {}, dtype::is_complex<T>()){}

        explicit VectorWriteManager(const v_t<T>* vec): VectorWriteManager(vec, vec->size()){}

        void* transfer(uint_t& nitem) override {
            nitem = nitem_next();
            const auto ptr = m_vec->data() + m_ncall_transfer * m_nelem_per_transfer;
            ++m_ncall_transfer;
            if (!nitem) return nullptr;
            return {reinterpret_cast<const void*>(ptr), nitem * m_item_size};
        }
    };

    /**
     * The chunk of data transacted in a single
     */
    struct Hyperslab {
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

        Hyperslab(hid_t parent_handle, str_t name, const uintv_t &item_dims, const uint_t &nitem,
                  bool writemode, hid_t h5type);

        Hyperslab(const IoManager& iom, str_t name, bool writemode);

        /**
         * sets the internal state to select the iitem-th item on this rank
         * @param iitem
         *  item index for selection
         */
        void select_hyperslab(uint_t iitem_begin, uint_t nitem);

        ~Hyperslab();
    };


#if 0
    struct WriteManager {
        /**
         * native type of the HDF5 dataset to which the data is to be written
         */
        const Type m_h5_type;
        /**
         * elemental layout of an item. an element is either a single value of type m_h5_type, or a pair (complex)
         */
        const uintv_t m_item_shape;
        /**
         * size in bytes of each item
         */
        const uint_t m_item_size;
        /**
         * total number of items (locally held)
         */
        const uint_t m_nitem;
        /**
         * native-typed layout of an item (takes a minor index of 2 into account if the data is complex valued)
         */
        const uintv_t m_native_item_shape;
        /**
         * native-typed layout of all locally stored items
         */
        const uintv_t m_native_shape;
        /**
         * maximum number of items to flush in each call to next()
         */
        const uint_t m_nitem_per_flush;
        /**
         * maximum number of (real or complex) elements per call to next()
         */
        const uint_t m_nelement_per_flush;
    private:
        uintv_t make_native_item_shape(bool complex) const {
            auto shape = m_item_shape;
            if (complex) shape.push_back(2);
            return shape;
        }
        uintv_t make_native_shape() const {
            auto shape = m_native_item_shape;
            shape.insert(shape.begin(), m_nitem);
            return shape;
        }
    protected:
        uint_t m_nflush = 0ul;

        uint_t nitem_next() const {
            const uint_t first_item = m_nflush * m_nitem_per_flush;
            if (first_item < m_nitem) return std::min(m_nitem_per_flush, m_nitem-first_item);
            return 0ul;
        }

    public:
        WriteManager(Type type, uintv_t item_shape, uint_t nitem, uint_t nitem_per_flush, bool complex):
            m_h5_type(type), m_item_shape(item_shape), m_nitem(nitem), m_native_item_shape(make_native_item_shape(complex)),
            m_native_shape(make_native_shape()), m_nitem_per_flush(nitem_per_flush),
            m_nelement_per_flush(m_nitem_per_flush * ) {}

        WriteManager(const WriteManager& other):
            WriteManager(other.m_h5_type, other.m_nprim_per_item, other.m_nitem, other.m_nitem_per_flush){}

        WriteManager& operator=(const WriteManager& other) = delete;

        struct Selection {
            const void* m_ptr;
            const uint_t m_size;
        };

        virtual Selection next() = 0;
    };



//    struct WriteBufferer: WriteManager {
//
//    };
#endif
#endif
}


#endif //M7_HDF5_BUFFERER_H
