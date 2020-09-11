//
// Created by RJA on 11/09/2020.
//

#ifndef M7_NUMERICARRAYFIELD_H
#define M7_NUMERICARRAYFIELD_H

#include "NumericArrayElement.h"

template<typename T>
class NumericArrayField : public NumericField<T> {
protected:
    const size_t m_size;

public:
    NumericArrayField(Table *table, size_t nelement, size_t size, const std::string &description = "") :
            NumericField<T>(table, nelement * size, description), m_size(size) {
    }

    using Field::m_nelement;
    NumericArrayElement<T> operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) {
        ASSERT(ielement < m_nelement);
        return NumericArrayElement<T>(this, m_size, irow, isegment, ielement);
    }

    bool is_complex() const override { return consts::is_complex<T>(); }

    std::vector<T> to_vector(size_t irow, size_t isegment, size_t ielement){
        return (*this)(irow, isegment, ielement).to_vector();
    }

};


#endif //M7_NUMERICARRAYFIELD_H
