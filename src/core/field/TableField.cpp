//
// Created by RJA on 27/10/2020.
//

#include "TableField.h"
#include "src/core/table/Table.h"


TableField::TableField(TableX *table, FieldData field_data,
             size_t nelement, std::string description) :
        m_table(table), m_data(field_data), m_description(description), m_nelement(nelement),
        m_size(nelement * m_data.m_element_size), m_offset(m_table->add_field(this)) {
    m_data.m_details["number of elements"] = std::to_string(m_nelement);
    m_data.m_details["field size (bytes)"] = std::to_string(m_size);
    m_data.m_details["offset from row (bytes)"] = std::to_string(m_offset);
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
