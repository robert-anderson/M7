//
// Created by RJA on 27/10/2020.
//

#include "NdFieldBase.h"
#include "Table.h"

//NdFieldBaseX::NdFieldBaseX(TableX *table, size_t nelement, size_t element_size,
//                           const std::type_info &type_info, std::string description) :
//        FieldBaseX(element_size, type_info),
//        m_table(table), m_description(description), m_nelement(nelement),
//        m_size(nelement * element_size), m_offset(m_table->add_field(this)) {}


NdFieldBaseX::NdFieldBaseX(TableX *table, FieldBaseX &&field,
                           size_t nelement, std::string description) :
        FieldBaseX(field), m_table(table), m_description(description), m_nelement(nelement),
        m_size(nelement * m_element_size), m_offset(m_table->add_field(this)) {
    m_details["number of elements"] = std::to_string(m_nelement);
    m_details["field size (bytes)"] = std::to_string(m_size);
    m_details["offset from row (bytes)"] = std::to_string(m_offset);
}

char *NdFieldBaseX::begin(const size_t &irow) const {
    return m_table->begin(irow) + m_offset;
}