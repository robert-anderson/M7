//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_VECINDSTABLE_H
#define M7_VECINDSTABLE_H

#include "MappedTable.h"

template<typename row_t>
struct VecIndsTable : MappedTable<row_t> {

    /**
     * vector of index vectors
     */
    v_t<uintv_t> m_inds;

    /**
     * won't compile unless the row defines a key_field_t;
     */
    typedef typename row_fields::Key<row_t>::type key_field_t;

    VecIndsTable(const row_t &row, MappedTableOptions opts) : MappedTable<row_t>(row, opts) {}

    VecIndsTable(const row_t &row) : VecIndsTable(row, MappedTableOptions()) {}

    VecIndsTable &operator=(const VecIndsTable &other) {
        *this = static_cast<const MappedTable<row_t>&>(other);
        return *this;
    }

    VecIndsTable(const VecIndsTable &other) : VecIndsTable(other.m_row, other.m_mapping_opts) {
        m_inds = other.m_inds;
    }

    bool operator==(const VecIndsTable &other) const {
        if (m_inds != other.m_inds) return false;
        if (!MappedTable<row_t>::operator==(other)) return false;
        return true;
    }

    const uintv_t* lookup(const key_field_t &key) const {
        auto &res = MappedTable<row_t>::lookup(key);
        if (!res) return nullptr;
        return m_inds.data()+res.index();
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

    const uintv_t* append(const key_field_t &key, uint_t ind) {
        auto ptr = m_inds.data() + lookup_or_insert(key).index();
        ptr->push_back(ind);
        return ptr;
    }

    const uintv_t* append(const row_t &src, uint_t ind) {
        auto ptr = m_inds.data() + lookup_or_insert(src).index();
        ptr->push_back(ind);
        return ptr;
    }

    void resize(uint_t nrec, double factor = -1.0) override {
        MappedTable<row_t>::resize(nrec, factor);
        m_inds.resize(this->nrecord());
    }
};

#endif //M7_VECINDSTABLE_H
