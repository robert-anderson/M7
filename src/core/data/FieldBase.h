//
// Created by rja on 02/10/2020.
//

#ifndef M7_FIELDBASE_H
#define M7_FIELDBASE_H

#include <cstddef>
#include <typeindex>
#include <string>
#include <src/core/util/utils.h>
#include <map>
#include "Table_NEW.h"

class FieldBase {
protected:
    void set_offsets();
    Table_NEW *m_table;
public:
    size_t m_offset = 0ul;
    size_t m_bit_offset = 0ul;
    const size_t m_element_size;
    const size_t m_nelement;
    const size_t m_size;
    const size_t m_dsize;
    // identifier for the stored type
    const std::type_info &m_type_info;
    const std::string m_description;

    FieldBase(Table_NEW *table, size_t element_size, size_t nelement,
              const std::type_info &type_info, std::string description="");

    bool is_same_type_as(const FieldBase &other) const;


    /*
     * pointer to beginning of a field
     */
    char *begin(const size_t &irow) const;

    /*
     * pointer to beginning of a byte-addressable View
     */
    char *begin(const size_t &irow, const size_t &ielement) const;

    size_t back_offset() const { return m_offset + m_size; }

    virtual std::string element_to_string(size_t irow, size_t ielement) const = 0;

    std::string to_string(size_t irow) const {
        std::string res;
        for (size_t i = 0ul; i < m_nelement; ++i) res += element_to_string(irow, i) + " ";
        return res;
    }

    virtual std::map<std::string, std::string> details() const {
        return
                {
                        {"offset (bytes)", std::to_string(m_offset)},
                        {"total size (bytes)", std::to_string(m_size)},
                        {"total size (words)", std::to_string(m_dsize)},
                        {"element size (bytes)", std::to_string(m_element_size)},
                        {"number of elements", std::to_string(m_nelement)}
                };
    }

    class View {
    protected:
        const FieldBase *m_field = nullptr;
        char *m_ptr;

        inline void init(const FieldBase *field, const size_t &irow, const size_t &ielement) {
            m_field = field;
            m_ptr = field->begin(irow, ielement);
        }

        inline void init(const FieldBase *field, char *ptr) {
            m_field = field;
            m_ptr = ptr;
        }

        /*
         * For instantiation from a Field
         */
        View(const FieldBase *field, const size_t &irow, const size_t &ielement) {
            init(field, irow, ielement);
            ASSERT(m_ptr);
        }

        inline defs::data_t *dptr() const {
            ASSERT(m_ptr);
            return (defs::data_t *) m_ptr;
        }

        /*
         * Element size as stored in Field
         */
        const size_t &size() const {
            return m_field->m_size;
        }

        /*
         * Element size in defs::data_t words as stored in Field
         */
        const size_t &dsize() const {
            return m_field->m_dsize;
        }

        /*
         * strictly for use with buffered views, since the view contents can't be initialised until
         * the Field's ctor has been called
         */
        View() : m_field(nullptr), m_ptr(nullptr) {}

        View(const View &other) : View() {
            m_field = other.m_field;
            m_ptr = other.m_ptr;
        }

        View &operator=(const View &other) {
            ASSERT(m_field);
            ASSERT(other.m_field->m_size == m_field->m_size);
            std::copy(other.m_ptr, other.m_ptr + other.m_field->m_size, m_ptr);
            return *this;
        }

        virtual std::string to_string() const = 0;
    };
};


#endif //M7_FIELDBASE_H
