//
// Created by RJA on 11/09/2020.
//

#include "PermanentField.h"


PermanentElement::PermanentElement(NumericField<char> *field, char *begin) : NumericElement(field, begin) {}

NumericElement<char> &PermanentElement::operator=(const char &v) {
    return NumericElement::operator=(v);
}

PermanentField::PermanentField(Table *table, size_t nmode, size_t nboson_cutoff, const std::string &description) :
        NumericField(table, nmode, description), m_nboson_cutoff(nboson_cutoff) {}
