//
// Created by RJA on 27/10/2020.
//

#include "TableField.h"
#include "src/core/table/Table.h"


TableField::TableField(Table *table, FieldData field_data,
             size_t nelement, std::string description) :
        m_table(table), m_data(field_data), m_description(description), m_nelement(nelement),
        m_size(nelement * m_data.m_element_size), m_offset(m_table->add_field(this)) {
    m_data.m_details["number of elements"] = std::to_string(m_nelement);
    m_data.m_details["field size (bytes)"] = std::to_string(m_size);
    m_data.m_details["offset from row (bytes)"] = std::to_string(m_offset);
}

TableField::TableField(const TableField &other):
        m_table(other.m_table->m_last_copied),
        m_data(other.m_data),
        m_description(other.m_description),
        m_nelement(other.m_nelement),
        m_size(other.m_size),
        m_offset(other.m_offset){
    auto offset = m_table->add_field(this);
    if (offset!=m_offset) mpi::stop_all("Offset should match that of copied field");
}

char *TableField::begin(const size_t &irow) const {
    return m_table->begin(irow) + m_offset;
}

char *TableField::raw_ptr(const size_t &irow, const size_t &ielement) const {
    return begin(irow) + ielement * m_data.m_element_size;
}

bool TableField::is_same_type_as(const TableField &other) const {
    return m_data.m_type_info == other.m_data.m_type_info;
}
