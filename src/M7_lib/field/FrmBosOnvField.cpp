//
// Created by rja on 11/08/2021.
//

#include "FrmBosOnvField.h"

FrmBosOnvField::FrmBosOnvField(Row *row, BasisData bd, std::string name) :
        base_t(m_frm, m_bos),
        m_frm(row, bd.m_frm, prefix("fermions", name)),
        m_bos(row, bd.m_bos, prefix("bosons", name)),
        m_decoded(*this){}

FrmBosOnvField::FrmBosOnvField(Row *row, size_t nsite, size_t nmode, size_t nboson_max, std::string name) :
        FrmBosOnvField(row, {FrmBasisData(nsite), {nmode, nboson_max}}){}

FrmBosOnvField::FrmBosOnvField(const FrmBosOnvField &other) :
    base_t(m_frm, m_bos), m_frm(other.m_frm), m_bos(other.m_bos), m_decoded(*this){}

FrmBosOnvField &FrmBosOnvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_frm = inds.first;
    m_bos = inds.second;
    return *this;
}
