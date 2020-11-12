//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/field/FermionOnvSpecifier.h"
#include "src/core/field/Views.h"

struct DecodedDeterminant {
    const size_t m_nbit;
    const size_t m_element_dsize;
    defs::det_work m_inds{};
    size_t m_nind = 0ul;

    explicit DecodedDeterminant(const FermionOnvSpecifier& spec):
    m_nbit(spec.m_nbit), m_element_dsize(spec.m_ndataword) {}

    explicit DecodedDeterminant(const views::FermionOnv &view):
    DecodedDeterminant(view.spec()){}

    virtual void update(const views::FermionOnv &det_elem) = 0;
};

struct OccupiedOrbitals : DecodedDeterminant {
    OccupiedOrbitals(const FermionOnvSpecifier& spec);
    OccupiedOrbitals(const views::FermionOnv &view);
    OccupiedOrbitals(const views::FermiBosOnv &view);
    void update(const views::FermionOnv &view) override;
};

struct VacantOrbitals : DecodedDeterminant {
    VacantOrbitals(const FermionOnvSpecifier& spec);
    VacantOrbitals(const views::FermionOnv &view);
    VacantOrbitals(const views::FermiBosOnv &view);
    void update(const views::FermionOnv &view) override;
};

#endif //M7_DECODEDDETERMINANT_H
