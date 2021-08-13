//
// Created by rja on 11/08/2021.
//

#ifndef M7_MAETABLE_H
#define M7_MAETABLE_H

#include "src/core/field/Fields.h"

struct MaeRow : public Row {
    field::MaeInds m_inds;
    field::Numbers<defs::wf_t, 1> m_values;

    field::MaeInds &key_field() {
        return m_inds;
    };

    MaeRow(size_t exsig, size_t nvalue):
        m_inds(this, exsig), m_values(this, {nvalue}, "values"){}
};

struct SpecMomsRow : public Row {
    field::SpecMomInds m_inds;
    field::Numbers<defs::wf_t, 1> m_values;

    field::SpecMomInds &key_field() {
        return m_inds;
    };

    SpecMomsRow(size_t exsig, size_t nvalue):
        m_inds(this, exsig), m_values(this, {nvalue}, "values"){}
};






#endif //M7_MAETABLE_H
