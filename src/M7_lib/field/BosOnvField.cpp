//
// Created by rja on 27/09/2021.
//

#include "BosOnvField.h"

#include <utility>

BosOnvField::BosOnvField(Row *row, const sys::bos::Basis& basis, std::string name) :
        base_t(row, {{basis.m_nmode}, {"boson mode occupations"}}, std::move(name)),
        m_basis(basis), m_decoded(*this) {
}

BosOnvField::BosOnvField(Row *row, const sys::Basis &basis, std::string name) : BosOnvField(row, basis.m_bos, name) {
    basis.require_pure_bos();
}

BosOnvField::BosOnvField(Row *row, const sys::Sector &hs, std::string name) : BosOnvField(row, hs.basis(), name){}

BosOnvField::BosOnvField(const BosOnvField &other) :
    base_t(other), m_hs(other.m_hs), m_decoded(*this) {}

BosOnvField &BosOnvField::operator=(const defs::inds &inds) {
    DEBUG_ASSERT_EQ(inds.size(), nelement(), "Vector is not the correct size");
    for (size_t i = 0ul; i < inds.size(); ++i) (*this)[i] = inds[i];
    return *this;
}

void BosOnvField::set_ops(const defs::inds &iops) {
    zero();
    for (auto& iop: iops) (*this)[iop]++;
}

size_t BosOnvField::nboson() const {
    return this->sum();
}