//
// Created by rja on 26/07/2021.
//

#include "BosOnvConnection.h"

BosOpPair::BosOpPair(size_t imode, size_t nop) : m_imode(imode), m_nop(nop){
    DEBUG_ASSERT_GT(m_nop, 0ul, "mode with no associated operators should not appear in operator product")
}

BosOps::BosOps(size_t nmode) {
    m_pairs.reserve(nmode);
}

const std::vector<BosOpPair> &BosOps::pairs() const {
    return m_pairs;
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

BosOnvConnection::BosOnvConnection(size_t nmode) : m_ann(nmode), m_cre(nmode){}

void BosOnvConnection::clear() {
    m_ann.clear();
    m_cre.clear();
}

size_t BosOnvConnection::size() const {
    return m_ann.size()+m_cre.size();
}

void BosOnvConnection::connect(const fields::BosOnv &src, const fields::BosOnv &dst) {
    DEBUG_ASSERT_EQ(src.nelement(), dst.nelement(), "src and dst ONVs are incompatible");
    clear();
    for (size_t imode=0ul; imode<src.nelement(); ++imode){
        if (src[imode] > dst[imode]) m_ann.add({imode, size_t(src[imode]-dst[imode])});
        else if (src[imode] < dst[imode]) m_cre.add({imode, size_t(dst[imode]-src[imode])});
    }
}

void BosOnvConnection::connect(const fields::BosOnv &src, const fields::BosOnv &dst, BosOps &com) {
    DEBUG_ASSERT_EQ(src.nelement(), dst.nelement(), "src and dst ONVs are incompatible");
    clear();
    for (size_t imode=0ul; imode<src.nelement(); ++imode){
        if (src[imode] > dst[imode]) m_ann.add({imode, size_t(src[imode]-dst[imode])});
        else if (src[imode] < dst[imode]) m_cre.add({imode, size_t(dst[imode]-src[imode])});
        else com.add({imode, src[imode]});
    }
}

void BosOnvConnection::apply(const fields::BosOnv &src, fields::BosOnv &dst) const {
    DEBUG_ASSERT_EQ(src.nelement(), dst.nelement(), "src and dst ONVs are incompatible");
    dst = src;
    for (auto& pair: m_ann.pairs()) dst[pair.m_imode]-=pair.m_nop;
    for (auto& pair: m_cre.pairs()) dst[pair.m_imode]+=pair.m_nop;
}

void BosOnvConnection::apply(const fields::BosOnv &src, BosOps &com) const {
    auto ann_iter = m_ann.pairs().cbegin();
    auto ann_end = m_ann.pairs().cend();
    for (size_t imode=0ul; imode<src.nelement(); ++imode) {
        if (ann_iter!=ann_end && ann_iter->m_imode==imode)
            com.add({imode, src[imode]-ann_iter->m_nop});
        else com.add({imode, src[imode]});
    }
}

void BosOnvConnection::apply(const fields::BosOnv &src, fields::BosOnv &dst, BosOps &com) const {
    apply(src, dst);
    apply(src, com);
}
