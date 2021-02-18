//
// Created by rja on 05/11/2020.
//

#ifndef M7_FERMIBOSONVENUMERATOR_H
#define M7_FERMIBOSONVENUMERATOR_H

#include "FermionOnvEnumerator.h"
#include "BosonOnvEnumerator.h"
#include "src/core/fieldz/FieldsZ.h"
#include "src/core/fieldz/BufferedFields.h"

class FermiBosOnvEnumerator : public Enumerator<fieldsz::Onv<1>> {
    FermionOnvEnumerator m_det_enum;
    BosonOnvEnumerator m_bonv_enum;
    buffered::FermionOnv m_fonv;
    buffered::BosonOnv m_bonv;
public:
    FermiBosOnvEnumerator(size_t nsite, size_t nelec, size_t nmode, size_t nboson_cutoff):
            m_det_enum(nsite, nelec), m_bonv_enum(nmode, nboson_cutoff), m_fonv(nsite), m_bonv(nmode){
        m_det_enum.next(m_fonv);
    }
    FermiBosOnvEnumerator(size_t nsite, size_t nelec, int spin, size_t nmode, size_t nboson_cutoff):
            m_det_enum(nsite, nelec, spin), m_bonv_enum(nmode, nboson_cutoff), m_fonv(nsite), m_bonv(nmode){
        m_det_enum.next_element(m_fonv);
    }

    bool next_element(fieldsz::Onv<1> &result) override {
        bool inner_allfound = !m_bonv_enum.next(m_bonv);
        if (inner_allfound) {
            m_bonv_enum.next(m_bonv);
            if (!m_det_enum.next(m_fonv)) return false;
        }
        result.m_fonv = m_fonv;
        result.m_bonv = m_bonv;
        return true;
    }
};


#endif //M7_FERMIBOSONVENUMERATOR_H
