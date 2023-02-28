//
// Created by rja on 26/01/23.
//

#ifndef M7_GLOBALACCUMULATION_H
#define M7_GLOBALACCUMULATION_H

#include "BufferedTable.h"

template<typename row_t>
class GlobalAccumulation {
    typedef typename row_fields::Key<row_t>::type key_field_t;
    typedef typename row_fields::Value<row_t>::type value_field_t;

    buffered::MappedTable<row_t> m_current;
    buffered::Table<row_t> m_deltas;
    buffered::Table<row_t> m_all_deltas;
    uint_t m_naccum = 0ul;

public:
    GlobalAccumulation(str_t name, const row_t& row):
        m_current(name+" current", row), m_deltas(name+" deltas", row), m_all_deltas(name+" gathered deltas", row){
        m_current.resize(1000ul);
        m_deltas.resize(100ul);
        m_all_deltas.resize(m_deltas.capacity()*mpi::nrank());
    }

    void add(const key_field_t& key, const value_field_t& value) {
        static_cast<Row&>(m_deltas.m_row).push_back_jump();
        row_fields::key(m_deltas.m_row) = key;
        row_fields::value(m_deltas.m_row) = value;
    }

    void update() {
        m_all_deltas.all_gatherv(m_deltas);
        m_deltas.clear();
        const auto& delta = m_all_deltas.m_row;
        const auto& key = row_fields::key(delta);
        const auto& value = row_fields::value(delta);
        for (delta.restart(); delta; ++delta){
            auto& row = m_current.lookup(key);
            if (!row) m_current.insert(key, row);
            row_fields::value(row) += value;
        }
        ++m_naccum;
    }

    const MappedTable<row_t>& current() const {
        return m_current;
    }

    uint_t naccum() const {
        return m_naccum;
    }

};


#endif //M7_GLOBALACCUMULATION_H
