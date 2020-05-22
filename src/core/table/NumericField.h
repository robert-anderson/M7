//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_NUMERICFIELD_H
#define SANDBOX2_NUMERICFIELD_H

#include "Field.h"
#include "NumericElement.h"
#include "src/core/util/defs.h"

template<typename T>
class NumericField : public Field {
public:
    NumericField(Table *table, size_t nelement = 1, const std::string &description = "") :
        Field(table, sizeof(T), nelement, typeid(T), description) {}

    NumericElement<T> operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) {
        ASSERT(ielement < m_nelement);
        return NumericElement<T>(this, begin(irow, isegment) + ielement * m_element_size);
    }

    bool is_complex() const override { return consts::is_complex<T>(); }

public:
    std::string to_string(size_t irow, size_t isegment, size_t ielement) override {
        return (*this)(irow, isegment, ielement).to_string();
    }

};

#endif //SANDBOX2_NUMERICFIELD_H
