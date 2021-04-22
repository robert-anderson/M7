//
// Created by rja on 01/04/2021.
//

#include "MevGroup.h"

FermionPromoter::FermionPromoter(size_t ncom, size_t nop_insert) :
        m_nop_insert(nop_insert),
        m_ncomb(integer_utils::combinatorial(ncom, nop_insert)),
        m_all_combs(nop_insert * m_ncomb) {
    if (!nop_insert) return;
    CombinationEnumerator ce(ncom, nop_insert);
    defs::inds comb(nop_insert);
    size_t icomb = ~0ul;
    while (ce.next(comb, icomb)) {
        for (size_t i = 0ul; i < nop_insert; ++i) {
            auto j = icomb * nop_insert + i;
            ASSERT(j < m_all_combs.size());
            m_all_combs[j] = comb[i];
        }
    }
}

const defs::mev_ind_t *FermionPromoter::begin(const size_t &icomb) const {
    ASSERT(icomb < m_ncomb);
    return m_all_combs.data() + icomb * m_nop_insert;
}

bool FermionPromoter::apply(const size_t &icomb, const conn::Antisym<> &conn, fields::FermionMevInds &inds) const {
    auto comb_begin = begin(icomb);
    inds.zero();
    size_t ann_passed = 0ul;
    size_t cre_passed = 0ul;
    for (size_t iins = 0ul; iins < m_nop_insert; ++iins) {
        auto ins = conn.com(comb_begin[iins]);
        while (ann_passed < conn.nann() && conn.ann(ann_passed) < ins) {
            inds.m_ann[ann_passed + iins] = conn.ann(ann_passed);
            ++ann_passed;
        }
        inds.m_ann[ann_passed + iins] = ins;

        while (cre_passed < conn.ncre() && conn.cre(cre_passed) < ins) {
            inds.m_cre[cre_passed + iins] = conn.cre(cre_passed);
            ++cre_passed;
        }
        inds.m_cre[cre_passed + iins] = ins;
    }
    auto phase = (ann_passed + cre_passed) & 1ul;
    phase = phase==conn.phase();

    // the rest of the promoted connection is the same as the connection
    // TODO: this with memcpy
    while (ann_passed < conn.nann()) {
        inds.m_ann[ann_passed + m_nop_insert] = conn.ann(ann_passed);
        ++ann_passed;
    }
    while (cre_passed < conn.ncre()) {
        inds.m_cre[cre_passed + m_nop_insert] = conn.cre(cre_passed);
        ++cre_passed;
    }
    return phase;
}
