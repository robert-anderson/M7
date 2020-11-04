//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/field/DeterminantSpecifier.h"
#include "src/core/field/Views.h"

struct DecodedDeterminant {
    const size_t m_nbit;
    const size_t m_element_dsize;
    defs::det_work m_inds{};
    size_t m_nind = 0ul;

    explicit DecodedDeterminant(const DeterminantSpecifier& spec):
    m_nbit(spec.m_nbit), m_element_dsize(spec.m_ndataword) {}

    explicit DecodedDeterminant(const views::Determinant &view):
    DecodedDeterminant(view.field()){}

    virtual void update(const views::Determinant &det_elem) = 0;
};

struct OccupiedOrbitals : DecodedDeterminant {
    OccupiedOrbitals(const DeterminantSpecifier& spec);
    OccupiedOrbitals(const views::Determinant &view);
    void update(const views::Determinant &view) override;
};

struct VacantOrbitals : DecodedDeterminant {
    VacantOrbitals(const DeterminantSpecifier& spec);
    VacantOrbitals(const views::Determinant &view);
    void update(const views::Determinant &view) override;
};

#endif //M7_DECODEDDETERMINANT_H
