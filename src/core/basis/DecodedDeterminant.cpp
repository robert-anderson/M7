//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

OccupiedOrbitals::OccupiedOrbitals(const DeterminantSpecifier &spec) : DecodedDeterminant(spec) {}

OccupiedOrbitals::OccupiedOrbitals(const views::Determinant &view) : DecodedDeterminant(view) {
    update(view);
}

void OccupiedOrbitals::update(const views::Determinant &view) {
    ASSERT(view.nbit() == m_nbit);
    ASSERT(view.ndataword() == m_element_dsize);
    m_nind = 0ul;
    for (size_t idataword = 0ul; idataword < view.ndataword(); ++idataword) {
        auto work = view.get_dataword(idataword);
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}

VacantOrbitals::VacantOrbitals(const DeterminantSpecifier &spec) : DecodedDeterminant(spec) {}

VacantOrbitals::VacantOrbitals(const views::Determinant &view) : DecodedDeterminant(view) {
    update(view);
}

void VacantOrbitals::update(const views::Determinant &view) {
    ASSERT(view.nbit() == m_nbit);
    ASSERT(view.ndataword() == m_element_dsize);
    m_nind = 0ul;
    for (size_t idataword = 0ul; idataword < view.ndataword(); ++idataword) {
        auto work = view.get_antidataword(idataword);
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}