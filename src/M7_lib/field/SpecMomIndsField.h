//
// Created by Robert J. Anderson on 12/08/2021.
//

#ifndef M7_SPECMOMINDSFIELD_H
#define M7_SPECMOMINDSFIELD_H

#include "RdmIndsField.h"

struct SpecMomIndsField : CompositeField<NumberField<rdm_ind_t>, NumberField<rdm_ind_t>> {
    typedef CompositeField<NumberField<rdm_ind_t>, NumberField<rdm_ind_t>> base_t;
    NumberField<rdm_ind_t> m_left, m_right;

    SpecMomIndsField(Row *row, str_t name = "indices");
};

#endif //M7_SPECMOMINDSFIELD_H
