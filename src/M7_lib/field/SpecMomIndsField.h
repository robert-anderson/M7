//
// Created by rja on 12/08/2021.
//

#ifndef M7_SPECMOMINDSFIELD_H
#define M7_SPECMOMINDSFIELD_H

#include "MaeIndsField.h"

struct SpecMomIndsField : CompositeField<MaeIndsField, MaeIndsField> {
    typedef CompositeField<MaeIndsField, MaeIndsField> base_t;
    const size_t m_exsig;
    MaeIndsField m_left, m_right;

    SpecMomIndsField(Row *row, size_t exsig, std::string name = "indices");

    SpecMomIndsField(const SpecMomIndsField &other);
};

#endif //M7_SPECMOMINDSFIELD_H