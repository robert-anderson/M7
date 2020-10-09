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
    const size_t m_dsize;
    // identifier for the stored type
    const std::type_info& m_type_info;
    FieldBase(Table* table, size_t element_size, size_t nelement, const std::type_info& type_info);

    bool is_same_type_as(const FieldBase& other) const;

    void set_offsets();

public:
    /*
     * pointer to beginning of a field
     */
    char* begin(const size_t& irow) const;
    /*
     * pointer to beginning of a byte-addressable View
     */
    char* begin(const size_t& irow, const size_t& ielement) const;
    const size_t& offset() const {return m_offset;}
    const size_t& size() const {return m_size;}
    size_t back_offset() const {return m_offset+m_size;}
    virtual std::string to_string(size_t irow) const = 0;

    class View {
    protected:
        const FieldBase *m_field = nullptr;
        char *m_ptr;

        inline void init(const FieldBase *field, const size_t &irow, const size_t &ielement) {
            m_field = field;
            m_ptr = field->begin(irow, ielement);
        }

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
        const size_t& size() const {
            return m_field->m_size;
        }
        /*
         * Element size in defs::data_t words as stored in Field
         */
        const size_t& dsize() const {
            return m_field->m_dsize;
        }
        /*
         * strictly for use with buffered views, since the view contents can't be initialised until
         * the Field's ctor has been called
         */
        View() : m_field(nullptr), m_ptr(nullptr){}
    };
};


#endif //M7_FIELDBASE_H
