//
// Created by rja on 05/11/2020.
//

#ifndef M7_FRMBOSONVENUMERATOR_H
#define M7_FRMBOSONVENUMERATOR_H

#include "FermionOnvEnumerator.h"
#include "BosonOnvEnumerator.h"
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedFields.h"

class FrmBosOnvEnumerator : public Enumerator<field::FrmBosOnv> {
    FermionOnvEnumerator m_fonv_enum;
    BosonOnvEnumerator m_bonv_enum;
    buffered::FrmOnv m_fonv;
    buffered::BosOnv m_bonv;
public:
    FrmBosOnvEnumerator(BasisDims bd, size_t nelec, size_t nboson_cutoff):
            m_fonv_enum(bd.m_nsite, nelec), m_bonv_enum(bd.m_nmode, nboson_cutoff),
            m_fonv(bd.frm()), m_bonv(bd.bos()){
        m_fonv_enum.next(m_fonv);
    }
    FrmBosOnvEnumerator(BasisDims bd, size_t nelec, int spin, size_t nboson_cutoff):
            m_fonv_enum(bd.m_nsite, nelec, spin),
            m_bonv_enum(bd.m_nmode, nboson_cutoff),
            m_fonv(bd.frm()), m_bonv(bd.bos()){
        m_fonv_enum.next_element(m_fonv);
    }

    bool next_element(field::FrmBosOnv &result) override {
        bool inner_allfound = !m_bonv_enum.next(m_bonv);
        if (inner_allfound) {
            m_bonv_enum.next(m_bonv);
            if (!m_fonv_enum.next(m_fonv)) return false;
        }
        result.m_frm = m_fonv;
        result.m_bos = m_bonv;
        return true;
    }
};


#endif //M7_FRMBOSONVENUMERATOR_H
