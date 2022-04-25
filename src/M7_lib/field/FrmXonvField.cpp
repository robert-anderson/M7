//
// Created by rja on 29/01/2022.
//

#include "FrmXonvField.h"

FrmXonvField::FrmXonvField(Row *row, const sys::frm::Basis& hs, std::string name) :
        CompositeField(m_ket, m_bra),
        m_ket(row, hs, prefix("ket", name)),
        m_bra(row, hs, prefix("bra", name)) {}

FrmXonvField::FrmXonvField(Row *row, const HilbertSpace& hs, std::string name) :
        FrmXonvField(row, hs.m_frm, name) {}

FrmXonvField::FrmXonvField(const FrmXonvField &other) :
        CompositeField(m_ket, m_bra), m_ket(other.m_ket), m_bra(other.m_bra) {}

FrmXonvField &FrmXonvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_ket = inds.first;
    m_bra = inds.second;
    return *this;
}