//
// Created by rja on 27/09/2021.
//

#include "BosOnvField.h"

BosOnvField::BosOnvField(Row *row, const BosBasisData& bd, std::string name) :
        base_t(row, {{bd.m_nmode}, {"boson mode occupations"}}, name),
        m_bd(bd), m_decoded(*this) {}

BosOnvField::BosOnvField(Row *row, const BasisData& bd, std::string name) : BosOnvField(row, bd.m_bos, name) {
    bd.require_pure_bos();
}

BosOnvField::BosOnvField(const BosOnvField &other) : base_t(other), m_bd(other.m_bd), m_decoded(*this) {}

BosOnvField &BosOnvField::operator=(const defs::inds &inds) {
    DEBUG_ASSERT_EQ(inds.size(), nelement(), "Vector is not the correct size");
    for (size_t i = 0ul; i < inds.size(); ++i) (*this)[i] = inds[i];
    return *this;
}

void BosOnvField::set_ops(const defs::inds &iops) {
    zero();
    for (auto& iop: iops) (*this)[iop]++;
}