//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

OccupiedOrbitals::OccupiedOrbitals(const FermionOnvSpecifier &spec) : DecodedDeterminant(spec) {}

OccupiedOrbitals::OccupiedOrbitals(const views::FermionOnv &view) : DecodedDeterminant(view) {
    update(view);
}

OccupiedOrbitals::OccupiedOrbitals(const views::FermiBosOnv &view): OccupiedOrbitals(view.m_fonv){}

void OccupiedOrbitals::update(const views::FermionOnv &view) {
    ASSERT(view.nbit() == m_nbit);
    ASSERT(view.ndataword() == m_element_dsize);
    m_nind = 0ul;
    for (size_t idataword = 0ul; idataword < view.ndataword(); ++idataword) {
        auto work = view.get_dataword(idataword);
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}

VacantOrbitals::VacantOrbitals(const FermionOnvSpecifier &spec) : DecodedDeterminant(spec) {}

VacantOrbitals::VacantOrbitals(const views::FermionOnv &view) : DecodedDeterminant(view) {
    update(view);
}

VacantOrbitals::VacantOrbitals(const views::FermiBosOnv &view) : VacantOrbitals(view.m_fonv){}

void VacantOrbitals::update(const views::FermionOnv &view) {
    ASSERT(view.nbit() == m_nbit);
    ASSERT(view.ndataword() == m_element_dsize);
    m_nind = 0ul;
    for (size_t idataword = 0ul; idataword < view.ndataword(); ++idataword) {
        auto work = view.get_antidataword(idataword);
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}