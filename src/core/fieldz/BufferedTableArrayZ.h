//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLEARRAYZ_H
#define M7_BUFFEREDTABLEARRAYZ_H

#include "TableZ.h"

template<typename row_t>
class BufferedTableArrayZ {
    typedef TableZ<row_t> table_t;
    Buffer m_buffer;
    std::vector<table_t> m_tables;

    const size_t &row_dsize() const {
        return static_cast<const TableBaseZ &>(m_tables[0]).m_row_dsize;
    }

    size_t window_dsize() const {
        return static_cast<const TableBaseZ &>(m_tables[0]).m_bw.dsize();
    }

    void update_nrow() {
        auto nrow = window_dsize() / row_dsize();
        for (size_t i = 0ul; i < ntable(); ++i) {
            static_cast<TableBaseZ &>(m_tables[i]).m_nrow = nrow;
        }
    }

public:

    const size_t &nrow_per_table() const {
        return static_cast<const TableBaseZ &>(m_tables[0]).m_nrow;
    }

    defs::data_t *dbegin() {
        return m_tables[0].dbegin();
    }

    const defs::data_t *dbegin() const {
        return m_tables[0].dbegin();
    }

    size_t buffer_dsize() const {
        return m_buffer.dsize();
    }

    size_t bw_dsize() const {
        return (*this)[0].bw_dsize();
    }

    BufferedTableArrayZ(std::string name, size_t ntable, row_t&& row):
            m_buffer(name, ntable) {
        m_tables.reserve(ntable);
        for (size_t itable = 0ul; itable < ntable; ++itable) {
            m_tables.emplace_back({row});
            static_cast<TableBaseZ &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }

    BufferedTableArrayZ(const BufferedTableArrayZ<table_t> &other) :
            m_buffer(other.m_buffer.m_name, other.ntable(), 0) {
        m_tables.reserve(other.ntable());
        for (size_t itable = 0ul; itable < other.ntable(); ++itable) {
            m_tables.emplace_back(other.m_tables[itable]);
            static_cast<TableBaseZ &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }

    size_t ntable() const { return m_tables.size(); }

    void resize(size_t nrow_per_table) {
        m_buffer.resize(ntable() * nrow_per_table * row_dsize());
        update_nrow();
    }

    void expand(size_t delta_nrow_per_table) {
        m_buffer.expand(ntable() * delta_nrow_per_table * row_dsize());
        update_nrow();
    }

    void expand(size_t delta_nrow_per_table, double expansion_factor) {
        m_buffer.expand(ntable() * delta_nrow_per_table * row_dsize(), expansion_factor);
        update_nrow();
    }

    table_t &operator[](const size_t &itable) {
        return m_tables[itable];
    }

    const table_t &operator[](const size_t &itable) const {
        return m_tables[itable];
    }

    defs::inds hwms() const {
        defs::inds res(ntable());
        for (size_t i = 0ul; i < ntable(); ++i) res[i] = (*this)[i].m_hwm;
        return res;
    }

    defs::inds displs() const {
        defs::inds res(ntable());
        for (size_t i = 0ul; i < ntable(); ++i) res[i] = bw_dsize() * i;
        return res;
    }

    void clear() {
        for (auto &table:m_tables) {
            static_cast<TableBaseZ &>(table).clear();
            static_cast<TableBaseZ &>(table).m_hwm = 0;
        }
    }

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }

    std::string to_string() {
        std::string tmp;
        for (size_t itable = 0ul; itable < ntable(); ++itable) {
            tmp+="Table array element " + std::to_string(itable) + ":\n";
            tmp+=static_cast<const TableZ<row_t> &>(m_tables[itable]).to_string()+"\n";
        }
        return tmp;
    }
};

#endif //M7_BUFFEREDTABLEARRAYZ_H