//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

OccupiedOrbitals::OccupiedOrbitals(const Field *field) : DecodedDeterminant(field) {}

OccupiedOrbitals::OccupiedOrbitals(const DeterminantElement &det_elem) : DecodedDeterminant(det_elem) {
    update(det_elem);
}

void OccupiedOrbitals::update(const DeterminantElement &det_elem) {
    assert(det_elem.nbit() == m_nbit);
    assert(det_elem.dsize() == m_element_dsize);
    m_nind = 0ul;
    size_t idataword = ~0ul;
    defs::data_t work;
    DeterminantElement::DatawordEnumerator enumerator(det_elem);
    while (enumerator.next(work, idataword)) {
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}

VacantOrbitals::VacantOrbitals(const Field *field) : DecodedDeterminant(field) {}

VacantOrbitals::VacantOrbitals(const DeterminantElement &det_elem) : DecodedDeterminant(det_elem) {
    update(det_elem);
}

void VacantOrbitals::update(const DeterminantElement &det_elem) {
    assert(det_elem.nbit() == m_nbit);
    assert(det_elem.dsize() == m_element_dsize);
    m_nind = 0ul;
    size_t idataword = ~0ul;
    defs::data_t work;
    DeterminantElement::AntiDatawordEnumerator enumerator(det_elem);
    while (enumerator.next(work, idataword)) {
        while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
    }
}