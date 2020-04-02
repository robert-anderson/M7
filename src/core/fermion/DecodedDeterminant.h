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
    size_t m_nind;
    const bool m_vacant;

    DecodedDeterminant(const DeterminantElement &det_elem, bool vacant=false):
    m_nbit(det_elem.nbit()), m_element_dsize(det_elem.dsize()),
    m_inds(m_nbit), m_vacant(vacant){}

    void update(const DeterminantElement &det_elem) {
        assert(det_elem.nbit()==m_nbit);
        assert(det_elem.dsize()==m_element_dsize);
        m_nind = 0ul;
        size_t idataword = ~0ul;
        defs::data_t work;
        DeterminantElement::DatawordEnumerator enumerator(det_elem, m_vacant);
        while (enumerator.next(work, idataword)){
            while (work) m_inds[m_nind++] = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
        }
    }
};

struct OccupiedOrbitals : DecodedDeterminant {
    OccupiedOrbitals(const DeterminantElement &det_elem);
};

struct VacantOrbitals : DecodedDeterminant {
    VacantOrbitals(const DeterminantElement &det_elem);
};

#endif //M7_DECODEDDETERMINANT_H
