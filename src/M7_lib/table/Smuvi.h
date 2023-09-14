//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_SMUVI_H
#define M7_SMUVI_H

#include "MappedTable.h"
#include "BufferedTable.h"

/**
 * Shared Memory Unordered map to Vectors of Indices(/Items)
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
