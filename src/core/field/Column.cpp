//
// Created by rja on 25/01/2021.
//

#include "Column.h"
#include "src/core/table/Table.h"

ColumnBase::ColumnBase(Table *table, ColumnData column_data,
                       size_t nelement, std::string description) :
        m_table(table), m_data(column_data), m_description(description), m_nelement(nelement),
        m_size(nelement * m_data.m_element_size), m_offset(m_table->add_column(this)) {
    m_data.m_details["number of elements"] = std::to_string(m_nelement);
    m_data.m_details["field size (bytes)"] = std::to_string(m_size);
    m_data.m_details["offset from row (bytes)"] = std::to_string(m_offset);
}

ColumnBase::ColumnBase(const ColumnBase &other):
        m_table(other.m_table->m_last_copied),
        m_data(other.m_data),
        m_description(other.m_description),
        m_nelement(other.m_nelement),
        m_size(other.m_size),
        m_offset(other.m_offset){
    auto offset = m_table->add_column(this);
    if (offset!=m_offset) mpi::stop_all("Offset should match that of copied field");
}

char *ColumnBase::raw_view(const size_t &irow, const size_t &ielement) const {
    return (char*)m_table->dbegin(irow) + m_offset + ielement * m_data.m_element_size;
}

bool ColumnBase::is_same_type_as(const ColumnBase &other) const {
    return m_data.m_type_info == other.m_data.m_type_info;
}
