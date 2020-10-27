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

    struct View {
        const FieldBaseX &m_field;
        char *m_ptr;

        View(const FieldBaseX& field, char* ptr): m_field(field), m_ptr(ptr){}

//        View(const FieldX &field, const size_t &irow, const size_t &iflat) :
//                m_field(field), m_ptr(field.raw_ptr(irow, iflat)) {}

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

    template<typename ...Args>
    std::pair<const char *, size_t> raw_view(const size_t &irow, Args...inds) const {
        return {raw_ptr(irow, inds...), m_element_size};
    }

    virtual std::string element_string(char* ptr) const = 0;

    virtual std::map<std::string, std::string> details() const {
        return {
            {"element size (bytes)", std::to_string(m_element_size)}
        };
    }
};


#endif //M7_FIELD_H
