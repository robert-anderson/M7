//
// Created by rja on 11/08/2021.
//

#ifndef M7_FRMBOSONVFIELD_H
#define M7_FRMBOSONVFIELD_H

#include "FrmOnvField.h"
#include "BosOnvField.h"
#include "CompositeField.h"
#include "src/core/caches/DecodedFrmBosOnv.h"

struct FrmBosOnvField : CompositeField<FrmOnvField, BosOnvField> {
    typedef CompositeField<FrmOnvField, BosOnvField> base_t;
    FrmOnvField m_frm;
    BosOnvField m_bos;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::FrmBosOnv m_decoded;
    FrmBosOnvField(Row* row, BasisData bd, std::string name="");

    FrmBosOnvField(const FrmBosOnvField& other);

    FrmBosOnvField& operator=(const FrmBosOnvField& other) {
        m_frm = other.m_frm;
        m_bos = other.m_bos;
        return *this;
    }

    FrmBosOnvField& operator=(const std::pair<defs::inds, defs::inds> &inds);

    const size_t &nsite() const {
        return m_frm.m_nsite;
    }

};



#endif //M7_FRMBOSONVFIELD_H
