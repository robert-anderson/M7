//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLEARRAY_H
#define M7_BUFFEREDTABLEARRAY_H

#include "MappedTable.h"

template<typename row_t, bool mapped=0>
class BufferedTableArray {
public:
    typedef typename std::conditional<mapped, MappedTable<row_t>, Table<row_t>>::type table_t;
private:
    Buffer m_buffer;
    std::vector<table_t> m_tables;

    const size_t &row_size() const {
        return static_cast<const TableBase &>(m_tables[0]).m_row_size;
    }

    size_t window_size() const {
        return static_cast<const TableBase &>(m_tables[0]).m_bw.size();
    }

    void update_nrow() {
        auto nrow = window_size() / row_size();
        for (size_t i = 0ul; i < ntable(); ++i) {
            static_cast<TableBase &>(m_tables[i]).m_nrow = nrow;
        }
    }

public:

    const size_t &nrow_per_table() const {
        return static_cast<const TableBase &>(m_tables[0]).m_nrow;
    }

    defs::buf_t *begin() {
        return m_tables[0].begin();
    }

    const defs::buf_t *begin() const {
        return m_tables[0].begin();
    }

    size_t buffer_size() const {
        return m_buffer.size();
    }

    size_t bw_size() const {
        return (*this)[0].bw_size();
    }

    BufferedTableArray(std::string name, size_t ntable, const table_t& table):
            m_buffer(name, ntable) {
        m_tables.reserve(ntable);
        for (size_t itable = 0ul; itable < ntable; ++itable) {
            m_tables.emplace_back(table);
            static_cast<TableBase &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }

    BufferedTableArray(const BufferedTableArray<table_t> &other) :
            m_buffer(other.m_buffer.m_name, other.ntable(), 0) {
        m_tables.reserve(other.ntable());
        for (size_t itable = 0ul; itable < other.ntable(); ++itable) {
            m_tables.emplace_back(other.m_tables[itable]);
            static_cast<TableBase &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }


    BufferedTableArray& operator=(const BufferedTableArray<row_t, mapped> &other) {
        m_buffer.resize(other.m_buffer.size());
        for (size_t itable = 0ul; itable < other.ntable(); ++itable) {
            auto& this_table = m_tables[itable];
            auto& other_table = other.m_tables[itable];
            this_table.clear();
            this_table.push_back(other_table.m_nrow);
            this_table.m_bw = other_table.m_bw;
        }
        return *this;
    }

    size_t ntable() const { return m_tables.size(); }

    void resize(size_t nrow_per_table, double factor=-1.0) {
        m_buffer.resize(ntable() * nrow_per_table * row_size(), factor);
        update_nrow();
    }

    void expand(size_t nrow_per_table) {
        resize(this->nrow_per_table()+nrow_per_table);
    }

    table_t &operator[](const size_t &itable) {
        DEBUG_ASSERT_LT(itable, ntable(), "Table array access OOB");
        return m_tables[itable];
    }

    const table_t &operator[](const size_t &itable) const {
        DEBUG_ASSERT_LT(itable, ntable(), "Table array access OOB");
        return m_tables[itable];
    }

    defs::inds hwms() const {
        defs::inds res(ntable());
        for (size_t i = 0ul; i < ntable(); ++i) res[i] = (*this)[i].m_hwm;
        return res;
    }

    defs::inds displs() const {
        defs::inds res(ntable());
        for (size_t i = 0ul; i < ntable(); ++i) res[i] = bw_size() * i;
        return res;
    }

    void clear() {
        for (auto &table:m_tables) {
            static_cast<TableBase &>(table).clear();
            static_cast<TableBase &>(table).m_hwm = 0;
        }
    }

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }

    std::string to_string() {
        std::string tmp;
        for (size_t itable = 0ul; itable < ntable(); ++itable) {
            tmp+="Table array element " + std::to_string(itable) + ":\n";
            tmp+= static_cast<const Table<row_t> &>(m_tables[itable]).to_string() + "\n";
        }
        return tmp;
    }
};

#endif //M7_BUFFEREDTABLEARRAY_H