//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_MAETABLE_H
#define M7_MAETABLE_H

#include <M7_lib/field/Fields.h>

struct RdmRow : public Row {
    field::RdmInds m_inds;
    field::Numbers<wf_t, 1> m_values;

    field::RdmInds &key_field();

    RdmRow(OpSig exsig, uint_t nvalue);
};

struct SpecMomsRow : public Row {
    field::SpecMomInds m_inds;
    field::Numbers<wf_t, 1> m_values;

    field::SpecMomInds &key_field();

    SpecMomsRow(uint_t nvalue);
};






#endif //M7_MAETABLE_H
