//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include "MappedTable.h"

template<typename row_t, typename table_impl_t>
class BufferedTable : public table_impl_t {
    static_assert(std::is_base_of<Table<row_t>, table_impl_t>::value, "Template args incompatible");
    Buffer m_buffer;
public:
    typedef table_impl_t table_t;
    using TableBase::m_bw;

    BufferedTable(str_t name, const table_t& table): table_t(table), m_buffer(name, 1ul) {
        TableBase::set_buffer(&m_buffer);
        ASSERT(static_cast<const Row&>(Table<row_t>::m_row).m_table);
    }

    BufferedTable(const table_t& table): BufferedTable("", table){}

    BufferedTable& operator=(const BufferedTable<row_t, table_t> &other) {
        table_t::operator=(other);
        return *this;
    }

    bool operator==(const BufferedTable<row_t, table_t> &other) const {
        return static_cast<const table_t&>(*this) == static_cast<const table_t&>(other);
    }

    BufferedTable(const BufferedTable<row_t, table_t> &other) : BufferedTable(other.m_buffer.m_name, other){
        *this = other;
        table_t::m_row.restart();
    }

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }

    double get_expansion_factor() const {
        return m_buffer.m_expansion_factor;
    }
};

namespace buffered {
    template <typename row_t>
    struct Table : BufferedTable<row_t, ::Table<row_t>> {
        Table(const row_t &row): BufferedTable<row_t, ::Table<row_t>>(Table<row_t>(row)){}
    };
    template <typename row_t>
    struct MappedTable : BufferedTable<row_t, ::MappedTable<row_t>> {
        MappedTable(const row_t &row, uint_t nbucket = 0ul, uint_t remap_nlookup = 0ul, double remap_ratio = 0.0):
            BufferedTable<row_t, ::Table<row_t>>(MappedTable<row_t>(row, nbucket, remap_nlookup, remap_ratio)){}
    };
}

#endif //M7_BUFFEREDTABLE_H
