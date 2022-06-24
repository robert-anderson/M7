//
// Created by Robert J. Anderson on 11/08/2021.
//

#include "FrmBosOnvField.h"

#include <utility>

FrmBosOnvField::FrmBosOnvField(
        Row *row, const sys::frm::Basis &frm_basis, const sys::bos::Basis &bos_basis, std::string name):
        base_t(m_frm, m_bos),
        m_frm(row, frm_basis, prefix("fermions", name)),
        m_bos(row, bos_basis, prefix("bosons", name)),
        m_decoded(*this) {
    DEBUG_ASSERT_TRUE(m_frm.m_basis == frm_basis, "fermion basis data incorrectly copied");
    DEBUG_ASSERT_TRUE(m_bos.m_basis == bos_basis, "boson basis data incorrectly copied");
}

FrmBosOnvField::FrmBosOnvField(Row *row, const sys::Basis &basis, std::string name) :
    FrmBosOnvField(row, basis.m_frm, basis.m_bos, std::move(name)){}

FrmBosOnvField::FrmBosOnvField(Row *row, const sys::Sector& sector, std::string name) :
        FrmBosOnvField(row, sector.basis(), std::move(name)){}

FrmBosOnvField::FrmBosOnvField(const FrmBosOnvField &other) :
    base_t(m_frm, m_bos), m_frm(other.m_frm), m_bos(other.m_bos), m_decoded(*this){}

FrmBosOnvField &FrmBosOnvField::operator=(const std::pair<uintv_t, uintv_t> &inds) {
    m_frm = inds.first;
    m_bos = inds.second;
    return *this;
}
