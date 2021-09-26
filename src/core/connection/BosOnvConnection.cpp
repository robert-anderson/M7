//
// Created by rja on 26/07/2021.
//

#include "BosOnvConnection.h"

BosOpPair::BosOpPair(size_t imode, size_t nop) : m_imode(imode), m_nop(nop){
    DEBUG_ASSERT_GT(m_nop, 0ul, "mode with no associated operators should not appear in operator product")
}

BosOps::BosOps(BasisDims bd){
    bd.require_pure_bos();
    m_pairs.reserve(bd.m_nmode);
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

const BosOpPair &BosOps::operator[](const size_t &ipair) const {
    DEBUG_ASSERT_LT(ipair, m_pairs.size(), "boson operator index OOB");
    return m_pairs[ipair];
}

BosOnvConnection::BosOnvConnection(BasisDims bd) :
    m_ann({0ul, bd.m_nmode}), m_cre({bd.m_nsite, 0ul}){}

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
