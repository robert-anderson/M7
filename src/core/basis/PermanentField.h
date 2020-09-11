//
// Created by RJA on 11/09/2020.
//

#ifndef M7_PERMANENTFIELD_H
#define M7_PERMANENTFIELD_H

#include "src/core/table/NumericField.h"

class PermanentField;

class PermanentElement : public NumericElement<char> {
public:
    PermanentElement(NumericField<char> *field, char *begin);

    NumericElement<char> &operator=(const char &v) override;
};

class PermanentField : public NumericField<char> {
    const size_t m_nboson_cutoff;
public:
    PermanentField(Table *table, size_t nelement, size_t nboson_cutoff, const std::string &description="");

    PermanentElement operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) {
        return PermanentElement(this, element_begin(irow, isegment, ielement));
    }

    const PermanentElement operator()(const size_t &irow, const size_t &isegment = 0, const size_t &ielement = 0) const;


};


#endif //M7_PERMANENTFIELD_H
