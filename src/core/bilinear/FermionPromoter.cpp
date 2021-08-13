//
// Created by rja on 10/08/2021.
//

#include <src/core/util/Foreach.h>
#include "FermionPromoter.h"


FermionPromoter::FermionPromoter(size_t ncom, size_t nop_insert) :
    m_nop_insert(nop_insert),
    m_ncomb(integer_utils::combinatorial(ncom, nop_insert)),
    m_all_combs(nop_insert * m_ncomb){
    if (!nop_insert) return;

    foreach::rtnd::Ordered<> foreach_comb(ncom, nop_insert);
    size_t icomb = 0ul;
    auto fn = [&]() {
        for (size_t i = 0ul; i < nop_insert; ++i) {
            auto j = icomb * nop_insert + i;
            ASSERT(j < m_all_combs.size());
            m_all_combs[j] = foreach_comb[i];
        }
        ++icomb;
    };
    foreach_comb(fn);
}

const defs::mev_ind_t *FermionPromoter::begin(const size_t &icomb) const {
    ASSERT(icomb < m_ncomb);
    return m_all_combs.data() + icomb * m_nop_insert;
}

bool FermionPromoter::apply(const size_t &icomb, const conn::FrmOnv &conn,
                            const FrmOps &com, MaeIndsPair &frm_inds) const {
    auto comb_begin = begin(icomb);
    frm_inds.zero();
    size_t ann_passed = 0ul;
    size_t cre_passed = 0ul;
    for (size_t iins = 0ul; iins < m_nop_insert; ++iins) {
        auto ins = com[comb_begin[iins]];
        while (ann_passed < conn.m_ann.size() && conn.m_ann[ann_passed] < ins) {
            frm_inds.m_ann[ann_passed + iins] = conn.m_ann[ann_passed];
            ++ann_passed;
        }
        frm_inds.m_ann[ann_passed + iins] = ins;

        while (cre_passed < conn.m_cre.size() && conn.m_cre[cre_passed] < ins) {
            frm_inds.m_cre[cre_passed + iins] = conn.m_cre[cre_passed];
            ++cre_passed;
        }
        frm_inds.m_cre[cre_passed + iins] = ins;
    }
    auto phase = (ann_passed + cre_passed) & 1ul;

    // the rest of the promoted connection is the same as the connection
    while (ann_passed < conn.m_ann.size()) {
        frm_inds.m_ann[ann_passed + m_nop_insert] = conn.m_ann[ann_passed];
        ++ann_passed;
    }
    while (cre_passed < conn.m_cre.size()) {
        frm_inds.m_cre[cre_passed + m_nop_insert] = conn.m_cre[cre_passed];
        ++cre_passed;
    }
    return phase;
}
