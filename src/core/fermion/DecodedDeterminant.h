//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/table/DeterminantField.h"

struct DecodedDeterminant {
    const size_t m_nbit;
    const size_t m_element_dsize;
    defs::inds m_inds;
    size_t m_nind = 0ul;

    explicit DecodedDeterminant(const Field* field):
    m_nbit(field->nbit()), m_element_dsize(field->element_dsize()),
    m_inds(m_nbit) {}

    explicit DecodedDeterminant(const DeterminantElement &det_elem):
    DecodedDeterminant(det_elem.field()){}

    virtual void update(const DeterminantElement &det_elem) = 0;
};

struct OccupiedOrbitals : DecodedDeterminant {
    OccupiedOrbitals(const Field* field);
    OccupiedOrbitals(const DeterminantElement &det_elem);
    void update(const DeterminantElement &det_elem) override;
};

struct VacantOrbitals : DecodedDeterminant {
    VacantOrbitals(const Field* field);
    VacantOrbitals(const DeterminantElement &det_elem);
    void update(const DeterminantElement &det_elem) override;
};

#endif //M7_DECODEDDETERMINANT_H