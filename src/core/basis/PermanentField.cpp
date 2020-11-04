//
// Created by RJA on 11/09/2020.
//

#if 0
#include "PermanentField.h"

PermanentField::PermanentField(Table *table, size_t nelement, size_t nmode, size_t nboson_cutoff, const std::string &description) :
        NumericArraySpecifier<uint8_t>(table, nelement, nmode, description), m_nboson_cutoff(nboson_cutoff) {}

PermanentElement::PermanentElement(PermanentField* field, const size_t &array_size, const size_t &irow, const size_t &isegment, const size_t &iarray) :
        NumericArrayElement<uint8_t>(field, array_size, irow, isegment, iarray) {}

#endif