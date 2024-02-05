//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include <utility>
#include "MappedTable.h"

template<typename row_t, typename table_impl_t>
class BufferedTable : public table_impl_t {
    static_assert(std::is_base_of<Table<row_t>, table_impl_t>::value, "Template args incompatible");
    Buffer m_buffer;
public:
    typedef table_impl_t table_t;
    using TableBase::m_bw;

    BufferedTable(str_t name, const table_t& table, Buffer::Permissions permissions): table_t(table),
        m_buffer(std::move(name), 1ul, permissions) {
        TableBase::set_buffer(&m_buffer);
        ASSERT(static_cast<const Row&>(Table<row_t>::m_row).m_table);
    }

    BufferedTable(const table_t& table, Buffer::Permissions permissions): BufferedTable("", table, permissions){}

    BufferedTable& operator=(const BufferedTable<row_t, table_t> &other) {
        table_t::operator=(other);
        return *this;
    }

    bool operator==(const BufferedTable<row_t, table_t> &other) const {
        return static_cast<const table_t&>(*this) == static_cast<const table_t&>(other);
    }

    BufferedTable(const BufferedTable<row_t, table_t> &other) :
        BufferedTable(other.m_buffer.m_name, other, other.m_buffer.m_permissions){
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
        Table(str_t name, const row_t &row, Buffer::Permissions permissions = {}):
            BufferedTable<row_t, ::Table<row_t>>(name, ::Table<row_t>(row), permissions){}
        Table(const row_t &row, Buffer::Permissions permissions = {}): Table("", row, permissions){}
    };
    template <typename row_t>
    struct MappedTable : BufferedTable<row_t, ::MappedTable<row_t>> {
        MappedTable(str_t name, const row_t &row, MappedTableOptions opts, Buffer::Permissions permissions = {}):
            BufferedTable<row_t, ::MappedTable<row_t>>(name, ::MappedTable<row_t>(row, opts), permissions){}
        MappedTable(const row_t &row, Buffer::Permissions permissions = {}):
            MappedTable("", row, permissions){}
        MappedTable(str_t name, const row_t &row, Buffer::Permissions permissions = {}):
            MappedTable(name, row, {}, permissions){}
        MappedTable(const row_t &row, MappedTableOptions opts, Buffer::Permissions permissions = {}):
            MappedTable("", row, opts, permissions){}
    };
}

#endif //M7_BUFFEREDTABLE_H
