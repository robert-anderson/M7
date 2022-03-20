//
// Created by rja on 29/01/2022.
//

#include "FrmXonvField.h"

FrmXonvField::FrmXonvField(Row *row, size_t nsite, std::string name) :
        CompositeField(m_ket, m_bra),
        m_ket(row, nsite, prefix("ket", name)),
        m_bra(row, nsite, prefix("bra", name)) {}

FrmXonvField::FrmXonvField(Row *row, BasisData bd, std::string name) :
        FrmXonvField(row, bd.m_nsite, name) {}

FrmXonvField::FrmXonvField(const FrmXonvField &other) :
        CompositeField(m_ket, m_bra), m_ket(other.m_ket), m_bra(other.m_bra) {}

FrmXonvField &FrmXonvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_ket = inds.first;
    m_bra = inds.second;
    return *this;
}