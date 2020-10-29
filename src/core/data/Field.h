//
// Created by rja on 21/10/2020.
//

#ifndef M7_FIELD_H
#define M7_FIELD_H

#include <cstddef>
#include <string>
#include <cstring>
#include <src/core/util/defs.h>
#include <map>

struct FieldBaseX {
    const size_t m_element_size;
    const std::type_info &m_type_info;
    std::map<std::string, std::string> m_details;

    struct View {
        const FieldBaseX &m_field;
        char *m_ptr;

        View(const FieldBaseX& field, char* ptr): m_field(field), m_ptr(ptr){}

        defs::data_t *dptr() const { return (defs::data_t *) m_ptr; }

        defs::data_t *dptr(const size_t &i) const {
            ASSERT(i * defs::nbyte_data < m_field.m_element_size);
            return ((defs::data_t *) m_ptr) + i;
        }

        virtual std::string to_string() const = 0;

    protected:

        View(const View &other) : m_field(other.m_field), m_ptr(other.m_ptr) {}

        View &operator=(const View &other) {
            ASSERT(m_field.m_element_size == other.m_field.m_element_size);
            if (&other != this)
                std::memcpy(m_ptr, other.m_ptr, m_field.m_element_size);
            return *this;
        }
    };

    FieldBaseX(size_t element_size, const std::type_info &type_info);

    bool is_same_type_as(const FieldBaseX& other) const;

    virtual std::string element_string(char* ptr) const = 0;

};


#endif //M7_FIELD_H
