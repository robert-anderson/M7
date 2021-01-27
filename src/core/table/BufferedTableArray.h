//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFEREDTABLEARRAY_H
#define M7_BUFFEREDTABLEARRAY_H


#include "Table.h"

template<typename table_t>
class BufferedTableArray {
    Buffer m_buffer;
    std::vector<table_t> m_tables;

    const size_t &row_dsize() const {
        return static_cast<const Table &>(m_tables[0]).m_row_dsize;
    }

    size_t window_dsize() const {
        return static_cast<const Table &>(m_tables[0]).m_bw.dsize();
    }

    const size_t &nrow_per_table() const {
        return static_cast<const Table &>(m_tables[0]).m_nrow;
    }

    void update_nrow() {
        auto nrow = window_dsize() / row_dsize();
        for (size_t i = 0ul; i < ntable(); ++i) {
            static_cast<Table &>(m_tables[i]).m_nrow = nrow;
        }
    }

public:

    defs::data_t *dbegin() {
        ASSERT(m_buffer.dbegin() == m_tables[0].dbegin());
        return m_buffer.dbegin();
    }

    const defs::data_t *dbegin() const {
        ASSERT(m_buffer.dbegin() == m_tables[0].dbegin());
        return m_buffer.dbegin();
    }

    size_t buffer_dsize() const {
        return m_buffer.dsize();
    }

    size_t bw_dsize() const {
        return (*this)[0].bw_dsize();
    }

    BufferedTableArray(std::string name, size_t ntable, const table_t& table): m_buffer(name, ntable, 0) {
        m_tables.reserve(ntable);
        for (size_t itable = 0ul; itable < ntable; ++itable) {
            m_tables.emplace_back(table);
            static_cast<Table &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }

    BufferedTableArray(const BufferedTableArray<table_t> &other) :
            m_buffer(other.m_buffer.m_name, other.ntable(), 0) {
        m_tables.reserve(other.ntable());
        for (size_t itable = 0ul; itable < other.ntable(); ++itable) {
            m_tables.emplace_back(other.m_tables[itable]);
            static_cast<Table &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }

    size_t ntable() const { return m_tables.size(); }

    void resize(size_t nrow_per_table) {
        m_buffer.resize(ntable() * nrow_per_table * row_dsize());
        update_nrow();
    }

    void expand(size_t delta_nrow) {
        m_buffer.expand(ntable() * delta_nrow * row_dsize());
        update_nrow();
    }

    void expand() {
        m_buffer.expand();
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
            static_cast<Table &>(table).clear();
            static_cast<Table &>(table).m_hwm = 0;
        }
    }

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }

    void print_contents() {
        for (size_t itable = 0ul; itable < ntable(); ++itable) {
            std::cout << "Table array element " << std::to_string(itable) << ":\n";
            static_cast<const Table &>(m_tables[itable]).print_contents();
            std::cout << std::endl;
        }
    }
};


#endif //M7_BUFFEREDTABLEARRAY_H
