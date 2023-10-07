//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_SMUVI_H
#define M7_SMUVI_H

#include "MappedTable.h"
#include "BufferedTable.h"

/**
 * Shared Memory Unordered map to Vectors of Indices(/Items) "SMUVI"
 * The SMUVI is a type of parallel hash map from a key domain to a variable-length vector of indices which has fewer
 * requirements than a parallelised generalization of a map from type T to a arbitrarily sized array of type U elements,
 * expressed using STL containers as std::unordered_map<T, std::vector<unsigned long>>.
 *
 * Primary among which is that the map is never simultaneously modified and read; it is at any given time in one of two
 * "insert" or "access"
 *
 * The required data structure need not support some operations commonly implemented in serial and concurrent hash maps:
 *  - Reassignment of existing elements. Newly inserted scalar values are appended to the existing list of values
 *    corresponding to the key
 *  - Deletion. No key or element of a value array may be removed after it has been added
 *  - Dynamic resizing of the "access" key array.
 *
 * The SMUVI behaves as a hash table which can be inserted to in parallel, then __collated__ with minimal
 * synchronisation overhead, to then be accessed for reading in constant time among all ranks in a communicator with a
 * common shared memory window.
 *
 * In the insert state, each rank initialises a private hash map and this is filled according to the workload of
 * insertions allotted to the rank. For access, shared-memory implementations of hash maps usually need to allow for the
 * possibility that a bin of keys can be written to or deleted from while the read operation is being carried out, and
 * as such a lock must be acquired before a read can take place, which results in sychronisation overhead. However, if
 * the SMUVI is in access state, each rank in a shared memory realm can lookup elements without the need to lock a
 * region of the buffer.
 *
 * The SMUVI state is transformed from insert to access in the collate operation. This proceeds by collecting the keys
 * from each insertion map into a local array which is then sorted. These sorted local keys are then gathered
 * contiguously into a shared memory buffer. The values are also transferred. The gathered keys are then globally sorted
 * via a K-way merge and duplicate keys are eliminated
 *
 * 1. each rank initializes a private MappedTable<row_t>, in which keys of the type row_t::key_t index the entries
 * 2. when a previously un-inserted
 *
 */

namespace smuvi {
    template<typename row_t>
    struct LocalTable : MappedTable<row_t> {
        /**
         * vector of index vectors
         */
        v_t<uintv_t> m_inds;

        /**
         * won't compile unless the row defines a key_field_t;
         */
        typedef typename row_fields::Key<row_t>::type key_field_t;

        LocalTable(const row_t &row, MappedTableOptions opts) : MappedTable<row_t>(row, opts) {}

        LocalTable(const row_t &row) : LocalTable(row, MappedTableOptions()) {}

        LocalTable &operator=(const LocalTable &other) {
            *this = static_cast<const MappedTable<row_t> &>(other);
            return *this;
        }

        LocalTable(const LocalTable &other) : LocalTable(other.m_row, other.m_mapping_opts) {
            m_inds = other.m_inds;
        }

        bool operator==(const LocalTable &other) const {
            if (m_inds != other.m_inds) return false;
            if (!MappedTable<row_t>::operator==(other)) return false;
            return true;
        }

        const uintv_t *inds(const row_t &row) const {
            return m_inds.data() + row.index();
        }

        const uintv_t *lookup(const key_field_t &key) const {
            auto &res = MappedTable<row_t>::lookup(key);
            if (!res) return nullptr;
            return inds(res);
        }

        void clear() override {
            MappedTable<row_t>::clear();
            m_inds.clear();
        }

        void erase(Lookup lookup) {
            m_inds.erase(m_inds.begin() + lookup.irow());
            MappedTable<row_t>::erase(lookup);
        }

        void erase(const key_field_t &key) {
            erase(lookup(key, MappedTable<row_t>::m_erase_row));
        }

        const uintv_t *append(const key_field_t &key, uint_t ind) {
            const auto i = MappedTable<row_t>::lookup_or_insert(key).index();
            auto ptr = m_inds.data() + i;
            ptr->push_back(ind);
            return ptr;
        }

        const uintv_t *append(const row_t &src, uint_t ind) {
            auto ptr = m_inds.data() + MappedTable<row_t>::lookup_or_insert(src).index();
            ptr->push_back(ind);
            return ptr;
        }

        void resize(uint_t nrec, double factor = -1.0) override {
            MappedTable<row_t>::resize(nrec, factor);
            m_inds.resize(this->m_bw.m_nrow);
        }
    };
}

namespace buffered {
    namespace smuvi {
        template<typename row_t>
        struct LocalTable : BufferedTable<row_t, ::smuvi::LocalTable<row_t>> {
            typedef BufferedTable<row_t, ::smuvi::LocalTable<row_t>> base_t;

            LocalTable(str_t name, const row_t &row, MappedTableOptions opts) :
                base_t(name, ::smuvi::LocalTable<row_t>(row, opts), false) {}

            LocalTable(const row_t &row): LocalTable("", row, {}) {}

            LocalTable(str_t name, const row_t &row) : LocalTable(name, row, {}) {}

            LocalTable(const row_t &row, MappedTableOptions opts): LocalTable("", row, opts) {}
        };
    }
}


//namespace smuvi {
//
//    template<typename row_t>
//    struct Smuvi {
//
//        void update(LocalTable<row_t>& local_table) {
//
//        }
//    };
//
//}

#endif //M7_SMUVI_H
