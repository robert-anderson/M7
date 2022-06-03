//
// Created by Robert J. Anderson on 27/09/2021.
//

#include "BosOnvField.h"
#include "M7_lib/connection/BosOnvConnection.h"

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
    base_t(other), m_basis(other.m_basis), m_decoded(*this) {}

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
size_t BosOnvField::occ_fac_square(const BosOnvConnection &conn) const {
    size_t fac = 1ul;
    for (auto& pair : conn.m_ann.pairs()) {
        for (size_t i=0ul; i<pair.m_nop; ++i) fac*= (*this)[pair.m_imode] - i;
    }
    for (auto& pair : conn.m_cre.pairs()) {
        for (size_t i=1ul; i<=pair.m_nop; ++i) fac*= (*this)[pair.m_imode] + i;
    }
    return fac;
}

double BosOnvField::occ_fac(const BosOnvConnection &conn) const {
    return std::sqrt(double(occ_fac_square(conn)));
}

size_t BosOnvField::occ_fac_square_ann(size_t occ, size_t nop) {
    size_t fac = 1ul;
    for (size_t i=0ul; i<nop; ++i){
        fac*=occ-i;
    }
    return fac;
}

size_t BosOnvField::occ_fac_square_cre(size_t occ, size_t nop) {
    size_t fac = 1ul;
    for (size_t i=1ul; i<=nop; ++i){
        fac*=occ+i;
    }
    return fac;
}

size_t BosOnvField::occ_fac_square_com(size_t occ, size_t nop) {
    size_t fac = 1ul;
    for (size_t i=0ul; i<nop; ++i){
        auto com_fac = occ-i;
        fac*=com_fac*com_fac;
    }
    return fac;
}

size_t BosOnvField::occ_fac_square(const BosOnvConnection &conn, const BosOps &com) const {
    size_t fac = 1;
    size_t icom = 0;
    auto ncom = com.pairs().size();
    /*
     * loop over bosonic annihilations
     */
    for (auto& pair : conn.m_ann.pairs()) {
        const size_t occ = (*this)[pair.m_imode];
        // if there are more annihilations than particles, the state is destroyed by the operator acting on it:
        if (pair.m_nop > occ) return 0;
        fac*= occ_fac_square_ann(occ, pair.m_nop);
        while (icom < ncom && com[icom].m_imode <= pair.m_imode){
            /*
             * there are common operator modes remaining and the current common mode index is not highter than the
             * current annihilation mode index
             */
            if (com[icom].m_imode == pair.m_imode) {
                // the current common operator mode is the same as the one just annihilated so deduct annihilated ops
                fac*= occ_fac_square_com(occ - pair.m_nop, com[icom].m_nop);
            }
            else {
                // the current common operator mode can't be in the annihilated array, since it's just been skipped over
                fac*= occ_fac_square_com((*this)[com[icom].m_imode], com[icom].m_nop);
            }
            ++icom;
        }
    }
    // only the annihilation part is affected by the common operators in normal ordering: do the creation ops as normal
    for (auto& pair : conn.m_cre.pairs()) fac*= occ_fac_square_cre((*this)[pair.m_imode], pair.m_nop);
    // do the rest of the common orbs that were not reached in the loop over annihilations
    for (;icom < ncom; ++icom) fac*= occ_fac_square_com((*this)[com[icom].m_imode], com[icom].m_nop);
    return fac;
}

double BosOnvField::occ_fac(const BosOnvConnection &conn, const BosOps &com) const {
    return std::sqrt(double(occ_fac_square(conn, com)));
}

size_t BosOnvField::occ_fac_square(const BosOps &com) const {
    size_t fac = 1ul;
    for (const auto& pair : com.pairs()) fac*= occ_fac_square_com((*this)[pair.m_imode], pair.m_nop);
    return fac;
}

double BosOnvField::occ_fac(const BosOps &com) const {
    return std::sqrt(double(occ_fac_square(com)));
}
