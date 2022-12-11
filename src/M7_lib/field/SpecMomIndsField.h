//
// Created by Robert J. Anderson on 12/08/2021.
//

#ifndef M7_SPECMOMINDSFIELD_H
#define M7_SPECMOMINDSFIELD_H

#include "MaeIndsField.h"

struct SpecMomIndsField : CompositeField<MaeIndsField, MaeIndsField> {
    typedef CompositeField<MaeIndsField, MaeIndsField> base_t;
    const uint_t m_exsig;
    MaeIndsField m_left, m_right;

    SpecMomIndsField(Row *row, OpSig exsig, str_t name = "indices");
};

#endif //M7_SPECMOMINDSFIELD_H
