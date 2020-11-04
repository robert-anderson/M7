//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/field/DeterminantSpecifier.h"

struct DecodedDeterminant {
    const size_t m_nbit;
    const size_t m_element_dsize;
    defs::det_work m_inds{};
    size_t m_nind = 0ul;

    explicit DecodedDeterminant(const DeterminantSpecifier& field):
    m_nbit(field.m_nbit), m_element_dsize(field.m_ndataword) {}

    explicit DecodedDeterminant(const DeterminantSpecifier::View &view):
    DecodedDeterminant(view.field()){}

    virtual void update(const DeterminantSpecifier::View &det_elem) = 0;
};

struct OccupiedOrbitals : DecodedDeterminant {
    OccupiedOrbitals(const DeterminantSpecifier& field);
    OccupiedOrbitals(const DeterminantSpecifier::View &view);
    void update(const DeterminantSpecifier::View &view) override;
};

struct VacantOrbitals : DecodedDeterminant {
    VacantOrbitals(const DeterminantSpecifier& field);
    VacantOrbitals(const DeterminantSpecifier::View &view);
    void update(const DeterminantSpecifier::View &view) override;
};

#endif //M7_DECODEDDETERMINANT_H
