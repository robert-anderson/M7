//
// Created by rja on 12/10/2020.
//

#ifndef M7_FORMATTEDFIELD_H
#define M7_FORMATTEDFIELD_H

#include "FieldBase.h"


struct SomeGetter {
    T& operator()(const FieldBase& field, const size_t& irow, const size_t& iflat){

    }
};




template<typename field_t>
class RowAccessor {
    char* m_row_begin = nullptr;
    field_t* m_field = nullptr;
};

template <size_t nind>
struct FormattedField : FieldBase {
    /*
    template<typename ...Args>
    FormattedField<nind>(Table_NEW *table, size_t element_size, const std::type_info &type_info, Args &&...shape) :
    FieldBase(table, element_size, nelement(shape...), type_info), m_format(std::forward<Args>(shape)...) {}
     */
    RowAccessor<char> operator[](const size_t& irow){

    }
};


#endif //M7_FORMATTEDFIELD_H
