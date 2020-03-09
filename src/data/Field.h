//
// Created by Robert John Anderson on 2020-03-09.
//

#ifndef M7_FIELD_H
#define M7_FIELD_H

#include <src/multidim/Indexer.h>
#include "FieldBase.h"

template<typename T, size_t nind = 1>
struct Field : FieldBase {
    const Indexer<nind> m_indexer;
    const size_t m_offset;

    template<typename ...Args>
    Field(FieldSet *field_set, Args... shape) :
        FieldBase(field_set, Indexer<nind>(shape...).nelement()),
        m_indexer(Indexer<nind>(shape...)), m_offset(field_set->add_field(this)) {}

private:
    T *flat_get(const size_t &irow, const size_t &flat) {
        assert(flat<m_length);
        return (T *) (m_field_set->m_buffer + irow * m_field_set->m_length + m_offset) + flat;
    }

public:
    template<typename ...Args>
    T *operator()(const size_t &irow, Args... inds) {
        return flat_get(m_indexer.get(inds...));
    }

    template<typename ...Args>
    void set(const T& value, const size_t &irow, Args... inds) {
        *(*this)(irow, inds...) = value;
    }

    template<typename ...Args>
    std::string to_string(size_t irow, Args... inds) {
        return std::to_string(*(*this)(irow, inds...));
    }
};


#endif //M7_FIELD_H
