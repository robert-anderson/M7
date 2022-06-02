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

defs::inds BosOps::get() const {
    defs::inds vec;
    vec.reserve(size());
    for (auto& pair: m_pairs) {
        for (size_t iop=0ul; iop<pair.m_nop; ++iop) vec.push_back(pair.m_imode);
    }
    DEBUG_ASSERT_EQ(vec.size(), size(), "incorrect number of operators in vector");
    return vec;
}

void BosOps::set(const defs::inds &imodes) {
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

void BosOps::set(size_t i) {
    clear();
    add({i, 1ul});
}

void BosOps::set(size_t i, size_t j) {
    clear();
    DEBUG_ASSERT_LE(i, j, "ops must be in ascending order");
    if (i==j) add({i, 2});
    else {
        add({i, 1});
        add({j, 1});
    }
}

void BosOps::set(size_t i, size_t j, size_t k) {
    clear();
    DEBUG_ASSERT_LE(i, j, "ops must be in ascending order");
    DEBUG_ASSERT_LE(j, k, "ops must be in ascending order");
    if (i==j) {
        if (j==k) {
            add({i, 3});
        }
        else {
            add({i, 2});
            add({k, 1});
        }
    }
    else {
        if (j==k) {
            add({i, 1});
            add({j, 2});
        }
        else {
            add({i, 1});
            add({j, 1});
            add({k, 1});
        }
    }
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