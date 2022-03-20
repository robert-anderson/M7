//
// Created by rja on 05/11/2020.
//

#ifndef M7_FRMBOSONVENUMERATOR_H
#define M7_FRMBOSONVENUMERATOR_H

#include <M7_lib/field/Fields.h>
#include <M7_lib/table/BufferedFields.h>

#include "FermionOnvEnumerator.h"
#include "BosonOnvEnumerator.h"

class FrmBosOnvEnumerator : public Enumerator<field::FrmBosOnv> {
    FermionOnvEnumerator m_fonv_enum;
    BosonOnvEnumerator m_bonv_enum;
    buffered::FrmOnv m_fonv;
    buffered::BosOnv m_bonv;
public:
    FrmBosOnvEnumerator(BasisData bd, size_t nelec, size_t nboson_cutoff):
            m_fonv_enum(bd.m_nsite, nelec), m_bonv_enum(bd.m_nmode, nboson_cutoff),
            m_fonv(bd.m_nsite), m_bonv(bd.m_nmode){
        m_fonv_enum.next(m_fonv);
    }
    FrmBosOnvEnumerator(BasisData bd, size_t nelec, int spin, size_t nboson_cutoff):
            m_fonv_enum(bd.m_nsite, nelec, spin),
            m_bonv_enum(bd.m_nmode, nboson_cutoff),
            m_fonv(bd.m_nsite), m_bonv(bd.m_nmode){
        m_fonv_enum.next_element(m_fonv);
    }

    bool next_element(field::FrmBosOnv &result) override {
        bool inner_allfound = !m_bonv_enum.next(m_bonv);
        if (inner_allfound) {
            m_bonv_enum.next(m_bonv);
            if (!m_fonv_enum.next(m_fonv)) return false;
        }
        result.m_frm = static_cast<const FrmOnvField&>(m_fonv);
        result.m_bos = static_cast<const BosOnvField&>(m_bonv);
        return true;
    }
};


#endif //M7_FRMBOSONVENUMERATOR_H
