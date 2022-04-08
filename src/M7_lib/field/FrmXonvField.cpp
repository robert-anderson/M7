//
// Created by rja on 29/01/2022.
//

#include "FrmXonvField.h"

FrmXonvField::FrmXonvField(Row *row, const FrmBasisData& frm, std::string name) :
        CompositeField(m_ket, m_bra),
        m_ket(row, frm.m_nsite, prefix("ket", name)),
        m_bra(row, frm.m_nsite, prefix("bra", name)) {}

FrmXonvField::FrmXonvField(Row *row, const BasisData& bd, std::string name) :
        FrmXonvField(row, bd.m_frm, name) {}

FrmXonvField::FrmXonvField(const FrmXonvField &other) :
        CompositeField(m_ket, m_bra), m_ket(other.m_ket), m_bra(other.m_bra) {}

FrmXonvField &FrmXonvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_ket = inds.first;
    m_bra = inds.second;
    return *this;
}