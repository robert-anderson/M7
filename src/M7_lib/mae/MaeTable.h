//
// Created by rja on 11/08/2021.
//

#ifndef M7_MAETABLE_H
#define M7_MAETABLE_H

#include <M7_lib/field/Fields.h>

struct MaeRow : public Row {
    field::MaeInds m_inds;
    field::Numbers<defs::wf_t, 1> m_values;

    field::MaeInds &key_field();

    MaeRow(size_t exsig, size_t nvalue);
};

struct SpecMomsRow : public Row {
    field::SpecMomInds m_inds;
    field::Numbers<defs::wf_t, 1> m_values;

    field::SpecMomInds &key_field();

    SpecMomsRow(size_t exsig, size_t nvalue);
};






#endif //M7_MAETABLE_H
