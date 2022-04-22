//
// Created by rja on 11/08/2021.
//

#ifndef M7_FRMBOSONVFIELD_H
#define M7_FRMBOSONVFIELD_H

#include <M7_lib/caches/DecodedFrmBosOnv.h>

#include "FrmOnvField.h"
#include "BosOnvField.h"
#include "CompositeField.h"

struct FrmBosOnvField : CompositeField<FrmOnvField, BosOnvField> {
    typedef CompositeField<FrmOnvField, BosOnvField> base_t;
    FrmOnvField m_frm;
    BosOnvField m_bos;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::FrmBosOnv m_decoded;
    FrmBosOnvField(Row* row, HilbertSpace hs, std::string name="");

    FrmBosOnvField(Row* row, const FrmHilbertSpace& frm_hs, const BosHilbertSpace& bos_hs, std::string name=""):
        FrmBosOnvField(row, HilbertSpace(frm_hs, bos_hs), name){}

    FrmBosOnvField(Row *row, size_t nsite, size_t nmode, size_t nboson_max=defs::max_bos_occ, std::string name = "");

    FrmBosOnvField(const FrmBosOnvField& other);

    FrmBosOnvField& operator=(const FrmBosOnvField& other) {
        m_frm = other.m_frm;
        m_bos = other.m_bos;
        return *this;
    }

    FrmBosOnvField& operator=(const std::pair<defs::inds, defs::inds> &inds);

};



#endif //M7_FRMBOSONVFIELD_H
