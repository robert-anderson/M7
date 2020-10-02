//
// Created by rja on 02/10/2020.
//

#ifndef M7_FIELDBASE_H
#define M7_FIELDBASE_H

#include <cstddef>
#include <typeindex>
#include <string>
#include <src/core/util/utils.h>
#include "Table.h"

class FieldBase {
protected:
    Table* m_table;
    size_t m_offset = 0ul;
    const size_t m_element_size;
    const size_t m_nelement;
    const size_t m_size;
    // identifier for the stored type
    const std::type_info& m_type_info;
    FieldBase(Table* table, size_t element_size, size_t nelement, const std::type_info& type_info);

    bool is_same_type_as(const FieldBase& other) const;

    void set_offsets();

public:
    char* begin(const size_t& irow) const;
    const size_t& offset() const {return m_offset;}
    const size_t& size() const {return m_size;}
    size_t back_offset() const {return m_offset+m_size;}
    virtual std::string to_string(size_t irow) const = 0;
};


#endif //M7_FIELDBASE_H
