//
// Created by rja on 05/11/2020.
//

#ifndef M7_FERMIBOSONVENUMERATOR_H
#define M7_FERMIBOSONVENUMERATOR_H

#include "FermionOnvEnumerator.h"
#include "BosonOnvEnumerator.h"
#include "src/core/field/Elements.h"

class FermiBosOnvEnumerator : public Enumerator<views::FbOnv> {
    FermionOnvEnumerator m_det_enum;
    BosonOnvEnumerator m_bonv_enum;
    elements::FermionOnv m_det;
    elements::BosonOnv m_bonv;
public:
    FermiBosOnvEnumerator(size_t nsite, size_t nelec, size_t nmode, size_t nboson_cutoff):
    m_det_enum(nsite, nelec), m_bonv_enum(nmode, nboson_cutoff), m_det(nsite), m_bonv(nmode){
        m_det_enum.next(m_det);
    }
    FermiBosOnvEnumerator(size_t nsite, size_t nelec, int spin, size_t nmode, size_t nboson_cutoff):
            m_det_enum(nsite, nelec, spin), m_bonv_enum(nmode, nboson_cutoff), m_det(nsite), m_bonv(nmode){
        m_det_enum.next_element(m_det);
    }

    bool next_element(views::FbOnv &result) override {
        bool inner_allfound = !m_bonv_enum.next(m_bonv);
        if (inner_allfound) {
            m_bonv_enum.next(m_bonv);
            if (!m_det_enum.next(m_det)) return false;
        }
        result.m_fonv = m_det;
        result.m_bonv = m_bonv;
        return true;
    }
};


#endif //M7_FERMIBOSONVENUMERATOR_H
