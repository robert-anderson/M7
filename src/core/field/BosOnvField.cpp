//
// Created by rja on 27/09/2021.
//

#include "BosOnvField.h"

BosOnvField::BosOnvField(Row *row, size_t nmode, std::string name) :
        base_t(row, {{nmode}, {"boson mode occupations"}}, name), m_nmode(m_nelement), m_decoded(*this) {}

BosOnvField::BosOnvField(Row *row, BasisData bd, std::string name) : BosOnvField(row, bd.m_nmode, name) {
    bd.require_pure_bos();
}

BosOnvField::BosOnvField(const BosOnvField &other) : base_t(other), m_nmode(m_nelement), m_decoded(*this) {}

BosOnvField &BosOnvField::operator=(const defs::inds &inds) {
    DEBUG_ASSERT_EQ(inds.size(), nelement(), "Vector is not the correct size");
    for (size_t i = 0ul; i < inds.size(); ++i) (*this)[i] = inds[i];
    return *this;
}