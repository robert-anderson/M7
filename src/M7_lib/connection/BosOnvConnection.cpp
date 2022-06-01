//
// Created by rja on 26/07/2021.
//

#include "BosOnvConnection.h"

BosOpPair::BosOpPair(size_t imode, size_t nop) : m_imode(imode), m_nop(nop){
    DEBUG_ASSERT_GT(m_nop, 0ul, "mode with no associated operators should not appear in operator product")
}

BosOps::BosOps(size_t nmode){
    m_pairs.reserve(nmode);
}

const std::vector<BosOpPair> &BosOps::pairs() const {
    return m_pairs;
}

defs::inds BosOps::to_vector() const {
    defs::inds vec;
    vec.reserve(size());
    for (auto& pair: m_pairs) {
        for (size_t iop=0ul; iop<pair.m_nop; ++iop) vec.push_back(pair.m_imode);
    }
    DEBUG_ASSERT_EQ(vec.size(), size(), "incorrect number of operators in vector");
    return vec;
}

void BosOps::from_vector(const defs::inds &imodes) {
    clear();
    if (imodes.empty()) return;
    size_t istart = 0;
    for(size_t i = 1; i<imodes.size(); ++i){
        if (imodes[istart]!=imodes[i]) {
            add({imodes[istart], i - istart});
            istart = i;
        }
    }
    add({imodes[istart], imodes.size()-istart});
}

void BosOps::clear() {
    m_pairs.clear();
    m_nop = 0ul;
}

size_t BosOps::size() const {
    return m_nop;
}

void BosOps::add(BosOpPair &&pair) {
    DEBUG_ASSERT_TRUE(m_pairs.empty() || pair.m_imode > m_pairs.back().m_imode,
                      "bosonic mode indices should be added in ascending order");
    DEBUG_ASSERT_LT(pair.m_imode, m_pairs.capacity(), "bosonic mode index is OOB");
    DEBUG_ASSERT_LT(m_pairs.size(), m_pairs.capacity(),
                    "should never have more distinct boson operators than modes");
    m_pairs.push_back(pair);
    m_nop += pair.m_nop;
}

void BosOps::set(BosOpPair &&pair) {
    clear();
    add(std::forward<BosOpPair>(pair));
}
void BosOps::set(const size_t &imode) {
    set({imode, 1ul});
}

void BosOps::set(const size_t &imode, const size_t &jmode) {
    clear();
    if (imode==jmode) add({imode, 2ul});
    else {
        add({imode, 1l});
        add({jmode, 1l});
    };
}

const BosOpPair &BosOps::operator[](const size_t &ipair) const {
    DEBUG_ASSERT_LT(ipair, m_pairs.size(), "boson operator index OOB");
    return m_pairs[ipair];
}

size_t BosOps::get_imode(size_t iop) const {
    DEBUG_ASSERT_LT(iop, m_nop, "Boson operator index OOB");
    for (const auto& pair: m_pairs) {
        if (iop<pair.m_nop) return pair.m_imode;
        iop-=pair.m_nop;
    }
    return ~0ul;
}
BosOnvConnection::BosOnvConnection(size_t nmode) : m_ann(nmode), m_cre(nmode){}

BosOnvConnection::BosOnvConnection(BasisData bd) : BosOnvConnection(bd.m_nmode){}

BosOnvConnection::BosOnvConnection(const BosOnvField &mbf) : BosOnvConnection(mbf.nelement()){}

void BosOnvConnection::clear() {
    m_ann.clear();
    m_cre.clear();
}

size_t BosOnvConnection::size() const {
    return m_ann.size()+m_cre.size();
}

void BosOnvConnection::connect(const BosOnvField &src, const BosOnvField &dst) {
    DEBUG_ASSERT_EQ(src.nelement(), dst.nelement(), "src and dst ONVs are incompatible");
    clear();
    for (size_t imode=0ul; imode<src.nelement(); ++imode){
        if (src[imode] > dst[imode]) m_ann.add({imode, size_t(src[imode]-dst[imode])});
        else if (src[imode] < dst[imode]) m_cre.add({imode, size_t(dst[imode]-src[imode])});
    }
}

void BosOnvConnection::connect(const BosOnvField &src, const BosOnvField &dst, BosOps &com) {
    DEBUG_ASSERT_EQ(src.nelement(), dst.nelement(), "src and dst ONVs are incompatible");
    clear();
    for (size_t imode=0ul; imode<src.nelement(); ++imode){
        if (src[imode] > dst[imode]) m_ann.add({imode, size_t(src[imode]-dst[imode])});
        else if (src[imode] < dst[imode]) m_cre.add({imode, size_t(dst[imode]-src[imode])});
        else com.add({imode, src[imode]});
    }
}

void BosOnvConnection::apply(const BosOnvField &src, BosOnvField &dst) const {
    DEBUG_ASSERT_EQ(src.nelement(), dst.nelement(), "src and dst ONVs are incompatible");
    dst = src;
    for (auto& pair: m_ann.pairs()) dst[pair.m_imode]-=pair.m_nop;
    for (auto& pair: m_cre.pairs()) dst[pair.m_imode]+=pair.m_nop;
}

void BosOnvConnection::apply(const BosOnvField &src, BosOps &com) const {
    auto ann_iter = m_ann.pairs().cbegin();
    auto ann_end = m_ann.pairs().cend();
    for (size_t imode=0ul; imode<src.nelement(); ++imode) {
        if (ann_iter!=ann_end && ann_iter->m_imode==imode)
            com.add({imode, src[imode]-ann_iter->m_nop});
        else com.add({imode, src[imode]});
    }
}

void BosOnvConnection::apply(const BosOnvField &src, BosOnvField &dst, BosOps &com) const {
    apply(src, dst);
    apply(src, com);
}

size_t BosOnvConnection::exsig() const {
    return exsig_utils::encode(0, 0, m_cre.size(), m_ann.size());
}

bool BosOnvConnection::respects_occ_range(const BosOnvField &src, size_t nboson_max) const {
    for (auto& pair: m_cre.pairs()) if (src[pair.m_imode]+pair.m_nop > nboson_max) return false;
    for (auto& pair: m_ann.pairs()) if (src[pair.m_imode] < pair.m_nop) return false;
    return true;
}

size_t BosOnvConnection::occ_fac_square(const BosOnvField &src) const {
    size_t fac = 1;
    for (auto& pair : m_ann.pairs()) {
        for (size_t i=0ul; i<pair.m_nop; ++i) fac*= src[pair.m_imode] - i;
    }
    for (auto& pair : m_cre.pairs()) {
        for (size_t i=1ul; i<=pair.m_nop; ++i) fac*= src[pair.m_imode] + i;
    }
    return fac;
}

double BosOnvConnection::occ_fac(const BosOnvField &src) const {
    return std::sqrt(double(occ_fac_square(src)));
}

size_t BosOnvConnection::occ_fac_square_ann(size_t occ, size_t nop) {
    size_t fac = 1;
    for (size_t i=0ul; i<nop; ++i){
        fac*=occ-i;
    }
    return fac;
}

size_t BosOnvConnection::occ_fac_square_cre(size_t occ, size_t nop) {
    size_t fac = 1;
    for (size_t i=1ul; i<=nop; ++i){
        fac*=occ+i;
    }
    return fac;
}

size_t BosOnvConnection::occ_fac_square_com(size_t occ, size_t nop) {
    size_t fac = 1;
    for (size_t i=0ul; i<nop; ++i){
        auto com_fac = occ-i;
        fac*=com_fac*com_fac;
    }
    return fac;
}

size_t BosOnvConnection::occ_fac_square(const BosOnvField &src, const BosOps &com) const {
    size_t fac = 1;
    size_t icom = 0;
    auto ncom = com.pairs().size();
    /*
     * loop over bosonic annihilations
     */
    for (auto& pair : m_ann.pairs()) {
        const size_t occ = src[pair.m_imode];
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
                fac*= occ_fac_square_com(src[com[icom].m_imode], com[icom].m_nop);
            }
            ++icom;
        }
    }
    // only the annihilation part is affected by the common operators in normal ordering: do the creation ops as normal
    for (auto& pair : m_cre.pairs()) fac*= occ_fac_square_cre(src[pair.m_imode], pair.m_nop);
    // do the rest of the common orbs that were not reached in the loop over annihilations
    for (;icom < ncom; ++icom) fac*= occ_fac_square_com(src[com[icom].m_imode], com[icom].m_nop);
    return fac;
}

double BosOnvConnection::occ_fac(const BosOnvField &src, const BosOps &com) const {
    return std::sqrt(double(occ_fac_square(src, com)));
}
