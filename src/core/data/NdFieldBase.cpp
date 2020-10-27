//
// Created by RJA on 27/10/2020.
//

#include "NdFieldBase.h"
#include "Table.h"

NdFieldBaseX::NdFieldBaseX(TableX *table, size_t nelement, size_t element_size,
                           const std::type_info &type_info, std::string description) :
        FieldBaseX(element_size, type_info),
        m_table(table), m_description(description), m_nelement(nelement),
        m_size(nelement * element_size), m_offset(m_table->add_field(this)) {}

char *NdFieldBaseX::begin(const size_t &irow) const {
    return m_table->begin(irow)+m_offset;
}

char *NdFieldBaseX::raw_ptr(const size_t &irow, const size_t &ielement) const {
    return begin(irow) + ielement * m_element_size;
}

std::map<std::string, std::string> NdFieldBaseX::details() const {
    auto map = FieldBaseX::details();
    map["number of elements"] = std::to_string(m_nelement);
    map["field size (bytes)"] = std::to_string(m_size);
    map["offset from row (bytes)"] = std::to_string(m_offset);
    return map;
}
