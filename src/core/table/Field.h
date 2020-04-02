//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_FIELD_H
#define SANDBOX2_FIELD_H

#include <cstddef>
#include <string>
#include <typeindex>

class Element;
class Table;

class Field {
protected:
    // the Table instance with which this Field is associated,
    // if this is nullptr, the field is disabled and takes up no space in a table
    Table *m_table;
    // the size of each element in bytes
    size_t m_element_size;
    // the size of each element in data_t (rounded up)
    size_t m_element_dsize;
    // the number of elements stored in the field
    const size_t m_nelement;
    // identifier for the stored type
    const std::type_index m_type_index;
    // byte offset from beginning of table row
    size_t m_offset;

    std::string m_description;

public:
    Field(Table *table, size_t element_size, size_t nelement, const std::type_info &type_info);

    Element element(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0);

    /*
    Element operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0){
        return element(irow, isegment, ielement);
    }
     */

    bool compatible_with(const Field &rhs) const;

    virtual size_t nbit();

    virtual size_t element_dsize() const;

    friend class Table;

    friend class Element;

protected:
    char *begin(const size_t &irow, const size_t &isegment = 0);

    char *element_begin(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0);

};


#endif //SANDBOX2_FIELD_H
