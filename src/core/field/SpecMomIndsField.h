//
// Created by rja on 12/08/2021.
//

#ifndef M7_SPECMOMINDSFIELD_H
#define M7_SPECMOMINDSFIELD_H

#include "MultiField.h"
#include "MaeIndsField.h"

struct SpecMomIndsField : MultiField<MaeIndsField, MaeIndsField> {
    typedef MultiField<MaeIndsField, MaeIndsField> base_t;
    const size_t m_exsig;
    const std::string m_name;
    MaeIndsField &m_left, m_right;

    SpecMomIndsField(Row *row, size_t exsig, std::string name = "indices");

    SpecMomIndsField(const SpecMomIndsField &other);
};


#if 0
#include "CompositeField.h"
#include "MaeIndsField.h"

struct SpecMomIndsField : CompositeField<MaeIndsField, MaeIndsField> {
    typedef CompositeField<MaeIndsField, MaeIndsField> base_t;
    const size_t m_exsig;
    MaeIndsField m_left, m_right;

    SpecMomIndsField(Row *row, size_t exsig, std::string name = "indices");

    SpecMomIndsField(const SpecMomIndsField &other);

    SpecMomIndsField& operator=(const SpecMomIndsField &other);
};
#endif

#endif //M7_SPECMOMINDSFIELD_H
