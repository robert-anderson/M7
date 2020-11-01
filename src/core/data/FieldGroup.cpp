//
// Created by RJA on 01/11/2020.
//

#include "FieldGroup.h"

size_t FieldGroup::add_field(const NdFieldBaseX *field) {
    m_fields.push_back(field);
    return 0ul;
}

size_t FieldGroup::nfield() const {
    return m_fields.size();
}
