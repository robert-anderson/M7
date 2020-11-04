//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

#if 0
OccupiedOrbitals::OccupiedOrbitals(const DeterminantSpecifier &field) : DecodedDeterminant(field) {}

OccupiedOrbitals::OccupiedOrbitals(const DeterminantSpecifier::View &view) : DecodedDeterminant(view) {
    update(view);
}

void OccupiedOrbitals::update(const DeterminantSpecifier::View &view) {
    ASSERT(view.nbit() == m_nbit);
    ASSERT(view.ndataword() == m_element_dsize);
    m_nind = 0ul;
    size_t idataword = ~0ul;
    defs::data_t work;
    DeterminantSpecifier::View::DatawordEnumerator enumerator(view);
    while (enumerator.next(work, idataword)) {
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}

VacantOrbitals::VacantOrbitals(const DeterminantSpecifier &field) : DecodedDeterminant(field) {}

VacantOrbitals::VacantOrbitals(const DeterminantSpecifier::View &view) : DecodedDeterminant(view) {
    update(view);
}

void VacantOrbitals::update(const DeterminantSpecifier::View &view) {
    ASSERT(view.nbit() == m_nbit);
    ASSERT(view.ndataword() == m_element_dsize);
    m_nind = 0ul;
    size_t idataword = ~0ul;
    defs::data_t work;
    DeterminantSpecifier::View::AntiDatawordEnumerator enumerator(view);
    while (enumerator.next(work, idataword)) {
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}
#endif