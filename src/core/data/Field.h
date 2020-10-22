//
// Created by rja on 21/10/2020.
//

#ifndef M7_FIELD_H
#define M7_FIELD_H

#include <cstddef>
#include <string>
#include <cstring>
#include <src/core/util/defs.h>

struct TableX;

struct FieldX {
    TableX *m_table;
    const size_t m_nelement;
    const size_t m_element_size;
    const size_t m_size;
    const size_t m_offset;
    const std::string m_description;


    struct View {
        const FieldX& m_field;
        char* m_ptr;
        View(const FieldX& field, char* ptr): m_field(field), m_ptr(ptr){}
        defs::data_t* dptr() const {return (defs::data_t*)m_ptr;}
        defs::data_t* dptr(const size_t& i) const {
            ASSERT(i*defs::nbyte_data < m_field.m_element_size);
            return ((defs::data_t*)m_ptr)+i;
        }
        virtual std::string to_string() const = 0;

    protected:

        View(const View& other): m_field(other.m_field), m_ptr(other.m_ptr){}
        View& operator=(const View& other){
            ASSERT(m_field.m_element_size == other.m_field.m_element_size);
            if (&other!=this)
                std::memcpy(m_ptr, other.m_ptr, m_field.m_element_size);
            return *this;
        }
    };


    FieldX(TableX *table, size_t nelement, size_t element_size, std::string description);

    char *begin(const size_t &irow) const;

    char *raw_ptr(const size_t &irow, const size_t &iflat) const;

    template<typename ...Args>
    std::pair<const char *, size_t> raw_view(const size_t &irow, Args...inds) const {
        return {raw_ptr(irow, inds...), m_element_size};
    }

    virtual std::string element_string(size_t irow, size_t ielement) const = 0;

    std::string to_string(size_t irow) const {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i)
            res+=" "+element_string(irow, i);
        return res;
    }
};


#endif //M7_FIELD_H
