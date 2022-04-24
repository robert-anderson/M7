//
// Created by rja on 11/08/2021.
//

#include "FrmBosOnvField.h"

FrmBosOnvField::FrmBosOnvField(Row *row, HilbertSpace hs, std::string name) :
        base_t(m_frm, m_bos),
        m_frm(row, hs.m_frm, prefix("fermions", name)),
        m_bos(row, hs.m_bos, prefix("bosons", name)),
        m_decoded(*this){
    DEBUG_ASSERT_TRUE(m_frm.m_hs==hs.m_frm, "fermion basis data incorrectly copied");
    DEBUG_ASSERT_TRUE(m_bos.m_hs==hs.m_bos, "boson basis data incorrectly copied");
}

FrmBosOnvField::FrmBosOnvField(const FrmBosOnvField &other) :
    base_t(m_frm, m_bos), m_frm(other.m_frm), m_bos(other.m_bos), m_decoded(*this){}

FrmBosOnvField &FrmBosOnvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_frm = inds.first;
    m_bos = inds.second;
    return *this;
}
