//
// Created by RJA on 11/09/2020.
//

#ifndef M7_PERMANENTFIELD_H
#define M7_PERMANENTFIELD_H

#include "src/core/table/NumericArrayField.h"

class PermanentField;

class PermanentElement : public NumericArrayElement<uint8_t> {
public:
    PermanentElement(PermanentField* field, const size_t &array_size, const size_t &irow, const size_t &isegment, const size_t &iarray);
};

class PermanentField : public NumericArrayField<uint8_t> {
    const size_t m_nboson_cutoff;
public:
    PermanentField(Table *table, size_t nelement, size_t nmode, size_t nboson_cutoff, const std::string &description="");

    PermanentElement operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) {
        return PermanentElement(this, m_size, irow, isegment, ielement);
    }

    const PermanentElement operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) const;

};


#endif //M7_PERMANENTFIELD_H
