//
// Created by rja on 11/08/2021.
//

#include "FrmBosOnvField.h"

FrmBosOnvField::FrmBosOnvField(Row *row, BasisData bd, std::string name) :
        base_t(m_frm, m_bos),
        m_frm(row, bd.m_nsite, prefix("fermions", name)),
        m_bos(row, bd.m_nmode, prefix("bosons", name)),
        m_decoded(*this){}

FrmBosOnvField::FrmBosOnvField(const FrmBosOnvField &other) :
    base_t(m_frm, m_bos), m_frm(other.m_frm), m_bos(other.m_bos), m_decoded(*this){}

FrmBosOnvField &FrmBosOnvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_frm = inds.first;
    m_bos = inds.second;
    return *this;
}
