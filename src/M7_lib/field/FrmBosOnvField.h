//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_FRMBOSONVFIELD_H
#define M7_FRMBOSONVFIELD_H

#include <M7_lib/caches/DecodedFrmBosOnv.h>

#include "FrmOnvField.h"
#include "BosOnvField.h"
#include "CompositeField.h"

struct FrmBosOnvField : CompositeField<FrmOnvField, BosOnvField> {
    typedef CompositeField<FrmOnvField, BosOnvField> base_t;
    const str_t m_name;
    FrmOnvField m_frm;
    BosOnvField m_bos;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::FrmBosOnv m_decoded;

    FrmBosOnvField(Row* row, const sys::frm::Basis& frm_basis, const sys::bos::Basis& bos_basis, str_t name="");
    /*
     * all ONVs implement the following ctor
     */
    FrmBosOnvField(Row* row, const sys::Basis& basis, str_t name="");
    /*
     * this particular MBF only needs the basis, but future MBF types might need the full sector information, and so
     * a common interface is realised by implementing a ctor of the following form in all MBFs
     */
    FrmBosOnvField(Row* row, const sys::Sector& sector, str_t name="");

    FrmBosOnvField(Row* row, const FrmBosOnvField& other);

    FrmBosOnvField(const FrmBosOnvField& other);

    FrmBosOnvField& operator=(const FrmBosOnvField& other) {
        m_frm = other.m_frm;
        m_bos = other.m_bos;
        return *this;
    }

    FrmBosOnvField& operator=(const std::pair<uintv_t, uintv_t>& inds);

};



#endif //M7_FRMBOSONVFIELD_H
