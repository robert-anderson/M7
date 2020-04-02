//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_NUMERICFIELD_H
#define SANDBOX2_NUMERICFIELD_H

#include <assert.h>
#include "Field.h"
#include "NumericElement.h"

template <typename T>
class NumericField : public Field{
public:
    NumericField(Table *table, size_t nelement=1):
    Field(table, sizeof(T), nelement, typeid(T)){}

    NumericElement<T> element(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) {
        assert(ielement<m_nelement);
        return NumericElement<T>(this, begin(irow, isegment)+ielement*m_element_size);
    }
};

#endif //SANDBOX2_NUMERICFIELD_H