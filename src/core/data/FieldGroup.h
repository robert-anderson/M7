//
// Created by RJA on 01/11/2020.
//

#ifndef M7_FIELDGROUP_H
#define M7_FIELDGROUP_H


#include "NdFieldBase.h"

struct FieldGroup {
    std::vector<const NdFieldBaseX *> m_fields;

    virtual size_t add_field(const NdFieldBaseX *field);

    size_t nfield() const;
};


#endif //M7_FIELDGROUP_H
