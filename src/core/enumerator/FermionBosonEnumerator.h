//
// Created by rja on 05/11/2020.
//

#ifndef M7_FERMIONBOSONENUMERATOR_H
#define M7_FERMIONBOSONENUMERATOR_H

#include "DeterminantEnumerator.h"
#include "BosonOnvEnumerator.h"
#include "src/core/field/Elements.h"

class FermionBosonEnumerator : public Enumerator<views::FermionBosonConfiguration> {
    DeterminantEnumerator m_det_enum;
    BosonOnvEnumerator m_bonv_enum;
    elements::Determinant m_det;
    elements::BosonOnv m_bonv;
public:
    FermionBosonEnumerator(size_t nsite, size_t nelec, size_t nmode, size_t occ_cutoff):
    m_det_enum(nsite, nelec), m_bonv_enum(nmode, occ_cutoff), m_det(nsite), m_bonv(nmode){
        m_det_enum.next(m_det);
    }
    FermionBosonEnumerator(size_t nsite, size_t nelec, int spin, size_t nmode, size_t occ_cutoff):
            m_det_enum(nsite, nelec, spin), m_bonv_enum(nmode, occ_cutoff), m_det(nsite), m_bonv(nmode){
        m_det_enum.next_element(m_det);
    }

    bool next_element(views::FermionBosonConfiguration &result) override {
        bool inner_allfound = !m_bonv_enum.next(m_bonv);
        if (inner_allfound) {
            m_bonv_enum.next(m_bonv);
            if (!m_det_enum.next(m_det)) return false;
        }
        result.m_det = m_det;
        result.m_perm = m_bonv;
        return true;
    }
};


#endif //M7_FERMIONBOSONENUMERATOR_H
