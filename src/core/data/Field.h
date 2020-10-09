//
// Created by rja on 02/10/2020.
//

#ifndef M7_FIELD_H
#define M7_FIELD_H

#include "FieldBase.h"
#include "src/core/nd/NdArrayFormat.h"

template<size_t nind>
class Field : public FieldBase {
protected:
    NdArrayFormat<nind> m_format;

    template<typename ...Args>
    size_t nelement(Args... args) {
        return NdArrayFormat<nind>(args...).nelement();
    }

    template<typename ...Args>
    Field(Table *table, size_t element_size, const std::type_info &type_info, Args &&...shape) :
            FieldBase(table, element_size, nelement(shape...), type_info), m_format(std::forward<Args>(shape)...) {}
};


#endif //M7_FIELD_H
