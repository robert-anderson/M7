//
// Created by rja on 05/11/2020.
//

#ifndef M7_FERMIBOSONVENUMERATOR_H
#define M7_FERMIBOSONVENUMERATOR_H

#include "FermionOnvEnumerator.h"
#include "BosonOnvEnumerator.h"
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedFields.h"

class FermiBosOnvEnumerator : public Enumerator<field::FrmBosOnv> {
    FermionOnvEnumerator m_det_enum;
    BosonOnvEnumerator m_bonv_enum;
    buffered::FrmOnv m_fonv;
    buffered::BosOnv m_bonv;
public:
    FermiBosOnvEnumerator(size_t nsite, size_t nelec, size_t nmode, size_t nboson_cutoff):
            m_det_enum(nsite, nelec), m_bonv_enum(nmode, nboson_cutoff), m_fonv(nsite), m_bonv(nmode){
        m_det_enum.next(m_fonv);
    }
    FermiBosOnvEnumerator(size_t nsite, size_t nelec, int spin, size_t nmode, size_t nboson_cutoff):
            m_det_enum(nsite, nelec, spin), m_bonv_enum(nmode, nboson_cutoff), m_fonv(nsite), m_bonv(nmode){
        m_det_enum.next_element(m_fonv);
    }

    bool next_element(field::FrmBosOnv &result) override {
        bool inner_allfound = !m_bonv_enum.next(m_bonv);
        if (inner_allfound) {
            m_bonv_enum.next(m_bonv);
            if (!m_det_enum.next(m_fonv)) return false;
        }
        result.m_frm = m_fonv;
        result.m_bos = m_bonv;
        return true;
    }
};


#endif //M7_FERMIBOSONVENUMERATOR_H
