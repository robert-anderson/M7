//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLEARRAY_H
#define M7_BUFFEREDTABLEARRAY_H

#include "MappedTable.h"

template<typename row_t, typename table_impl_t>
class BufferedTableArray {
    static_assert(std::is_base_of<Table<row_t>, table_impl_t>::value, "Template args incompatible");
    Buffer m_buffer;
    v_t<table_impl_t> m_tables;

    uint_t row_size() const {
        return static_cast<const TableBase&>(m_tables[0]).slot_size();
    }

    uint_t window_size() const {
        return static_cast<const TableBase &>(m_tables[0]).m_bw.m_size;
    }

public:
    typedef table_impl_t table_t;

    uint_t nrow_per_table() const {
        return static_cast<const TableBase&>(m_tables[0]).nslot();
    }

    buf_t *begin() {
        return m_tables[0].begin();
    }

    const buf_t *begin() const {
        return m_tables[0].begin();
    }

    uint_t buffer_size() const {
        return m_buffer.size();
    }

    uint_t bw_size() const {
        return (*this)[0].bw_size();
    }

    BufferedTableArray(str_t name, uint_t ntable, const row_t& row): m_buffer(name, ntable) {
        m_tables.reserve(ntable);
        for (uint_t itable = 0ul; itable < ntable; ++itable) {
            m_tables.emplace_back(table_t(row));
            static_cast<TableBase &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }

    BufferedTableArray(const BufferedTableArray<row_t, table_t> &other) :
            m_buffer(other.m_buffer.m_name, other.ntable()) {
        m_tables.reserve(other.ntable());
        for (uint_t itable = 0ul; itable < other.ntable(); ++itable) {
            m_tables.emplace_back(other.m_tables[itable]);
            static_cast<TableBase &>(m_tables.back()).set_buffer(&m_buffer);
        }
    }


    BufferedTableArray& operator=(const BufferedTableArray<row_t, table_t> &other) {
        m_buffer.resize(other.m_buffer.size());
        for (uint_t itable = 0ul; itable < other.ntable(); ++itable) {
            auto& this_table = m_tables[itable];
            auto& other_table = other.m_tables[itable];
            this_table.clear();
            this_table.push_back(other_table.m_nslot);
            this_table.m_bw = other_table.m_bw;
        }
        return *this;
    }

    uint_t ntable() const { return m_tables.size(); }

    void resize(uint_t nrow_per_table, double factor=-1.0) {
        m_buffer.resize(ntable() * nrow_per_table * row_size(), factor);
    }

    void expand(uint_t nrow_per_table) {
        resize(this->nrow_per_table()+nrow_per_table);
    }

    table_t &operator[](const uint_t &itable) {
        DEBUG_ASSERT_LT(itable, ntable(), "Table array access OOB");
        return m_tables[itable];
    }

    const table_t &operator[](const uint_t &itable) const {
        DEBUG_ASSERT_LT(itable, ntable(), "Table array access OOB");
        return m_tables[itable];
    }

    uintv_t hwms() const {
        uintv_t res(ntable());
        for (uint_t i = 0ul; i < ntable(); ++i) res[i] = (*this)[i].m_hwm;
        return res;
    }

    uintv_t displs() const {
        uintv_t res(ntable());
        for (uint_t i = 0ul; i < ntable(); ++i) res[i] = bw_size() * i;
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

    str_t to_string() {
        str_t tmp;
        for (uint_t itable = 0ul; itable < ntable(); ++itable) {
            tmp+="Table array element " + std::to_string(itable) + ":\n";
            tmp+= static_cast<const Table<row_t> &>(m_tables[itable]).to_string() + "\n";
        }
        return tmp;
    }
//
//    str_t name() const {
//        return m_buffer.m_name;
//    }
};


namespace buffered {
    template <typename row_t> using Tables = BufferedTableArray<row_t, ::Table<row_t>>;
    template <typename row_t> using MappedTables = BufferedTableArray<row_t, ::MappedTable<row_t>>;
}

#endif //M7_BUFFEREDTABLEARRAY_H