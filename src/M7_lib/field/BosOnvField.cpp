//
// Created by Robert J. Anderson on 27/09/2021.
//

#include "BosOnvField.h"
#include "M7_lib/connection/BosOnvConnection.h"

#include <utility>

BosOnvField::BosOnvField(Row *row, const sys::bos::Basis& basis, str_t name) :
        base_t(row, {{basis.m_nmode}, {"boson mode occupations"}}, std::move(name), true),
        m_basis(basis), m_decoded(*this) {
}

BosOnvField::BosOnvField(Row *row, const sys::Basis &basis, str_t name) : BosOnvField(row, basis.m_bos, name) {
    basis.require_pure_bos();
}

BosOnvField::BosOnvField(Row *row, const sys::bos::Sector &sector, str_t name) :
    BosOnvField(row, sector.m_basis, name){}

BosOnvField::BosOnvField(Row *row, const sys::Sector &sector, str_t name) :
    BosOnvField(row, sector.m_bos, name){}

BosOnvField::BosOnvField(Row* row, const BosOnvField& other) : BosOnvField(row, other.m_basis, other.m_name){}

BosOnvField::BosOnvField(const BosOnvField &other) :
    base_t(other), m_basis(other.m_basis), m_decoded(*this) {}

BosOnvField &BosOnvField::operator=(const uintv_t &inds) {
    DEBUG_ASSERT_EQ(inds.size(), nelement(), "Vector is not the correct size");
    for (uint_t i = 0ul; i < inds.size(); ++i) (*this)[i] = inds[i];
    return *this;
}

void BosOnvField::set_ops(const uintv_t &iops) {
    zero();
    for (auto& iop: iops) (*this)[iop]++;
}

uint_t BosOnvField::nboson() const {
    return this->sum();
}
uint_t BosOnvField::occ_fac_square(const BosOnvConnection &conn) const {
    uint_t fac = 1ul;
    for (auto& pair : conn.m_ann.pairs()) {
        for (uint_t i=0ul; i<pair.m_nop; ++i) fac*= (*this)[pair.m_imode] - i;
    }
    for (auto& pair : conn.m_cre.pairs()) {
        for (uint_t i=1ul; i<=pair.m_nop; ++i) fac*= (*this)[pair.m_imode] + i;
    }
    return fac;
}

double BosOnvField::occ_fac(const BosOnvConnection &conn) const {
    return std::sqrt(double(occ_fac_square(conn)));
}

uint_t BosOnvField::occ_fac_square_ann(uint_t occ, uint_t nop) {
    uint_t fac = 1ul;
    for (uint_t i=0ul; i<nop; ++i){
        fac*=occ-i;
    }
    return fac;
}

uint_t BosOnvField::occ_fac_square_cre(uint_t occ, uint_t nop) {
    uint_t fac = 1ul;
    for (uint_t i=1ul; i<=nop; ++i){
        fac*=occ+i;
    }
    return fac;
}

uint_t BosOnvField::occ_fac_square_com(uint_t occ, uint_t nop) {
    uint_t fac = 1ul;
    for (uint_t i=0ul; i<nop; ++i){
        auto com_fac = occ-i;
        fac*=com_fac*com_fac;
    }
    return fac;
}

uint_t BosOnvField::occ_fac_square(const BosOnvConnection &conn, const BosOps &com) const {
    uint_t fac = 1;
    uint_t icom = 0;
    auto ncom = com.pairs().size();
    /*
     * loop over bosonic annihilations
     */
    for (auto& pair : conn.m_ann.pairs()) {
        const uint_t occ = (*this)[pair.m_imode];
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

uint_t BosOnvField::occ_fac_square(const BosOps &com) const {
    uint_t fac = 1ul;
    for (const auto& pair : com.pairs()) fac*= occ_fac_square_com((*this)[pair.m_imode], pair.m_nop);
    return fac;
}

double BosOnvField::occ_fac(const BosOps &com) const {
    return std::sqrt(double(occ_fac_square(com)));
}

uint_t BosOnvField::occ_npair() const {
    uint_t n = 0ul;
    for (uint_t imode = 0ul; imode < m_basis.m_nmode; ++imode) {
        const uint_t occ = (*this)[imode];
        n += occ*(occ-1);
    }
    return n;
}