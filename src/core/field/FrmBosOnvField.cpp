//
// Created by rja on 11/08/2021.
//

#include "FrmBosOnvField.h"

FrmBosOnvField::FrmBosOnvField(Row *row, size_t nsite, std::string name) :
        MultiField<FrmOnvField, BosOnvField>(row,
                                             {nullptr, nsite, name.empty() ? "" : name + " (fermion)"},
                                             {nullptr, nsite, name.empty() ? "" : name + " (boson)"}),
        m_name(name), m_frm(get<0>()), m_bos(get<1>()) {
}

FrmBosOnvField::FrmBosOnvField(const FrmBosOnvField &other) : FrmBosOnvField(other.m_frm.row_of_copy(), other.m_frm.m_nsite,
                                                                             other.m_name) {}

FrmBosOnvField &FrmBosOnvField::operator=(const FrmBosOnvField &other) {
    m_frm = other.m_frm;
    m_bos = other.m_bos;
    return *this;
}

FrmBosOnvField &FrmBosOnvField::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    m_frm = inds.first;
    m_bos = inds.second;
    return *this;
}

const size_t &FrmBosOnvField::nsite() const {
    return m_frm.m_nsite;
}