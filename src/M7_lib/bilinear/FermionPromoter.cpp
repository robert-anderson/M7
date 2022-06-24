//
// Created by Robert J. Anderson on 10/08/2021.
//

#include <M7_lib/foreach/BasicForeach.h>
#include "FermionPromoter.h"


FermionPromoter::FermionPromoter(uint_t ncom, uint_t nop_insert) :
    m_nop_insert(nop_insert),
    m_ncomb(integer::combinatorial(ncom, nop_insert)),
    m_all_combs(nop_insert * m_ncomb){
    if (!nop_insert) return;

    basic_foreach::rtnd::Ordered<> foreach_comb(ncom, nop_insert);
    uint_t icomb = 0ul;
    auto fn = [&](const defs::uintv_t& inds) {
        for (uint_t i = 0ul; i < nop_insert; ++i) {
            auto j = icomb * nop_insert + i;
            ASSERT(j < m_all_combs.size());
            m_all_combs[j] = inds[i];
        }
        ++icomb;
    };
    foreach_comb.loop(fn);
}

const defs::mev_ind_t *FermionPromoter::begin(const uint_t &icomb) const {
    ASSERT(icomb < m_ncomb);
    return m_all_combs.data() + icomb * m_nop_insert;
}

bool FermionPromoter::apply(const uint_t &icomb, const conn::FrmOnv &conn,
                            const FrmOps &com, MaeIndsPair &frm_inds) const {
    auto comb_begin = begin(icomb);
    frm_inds.zero();
    uint_t ann_passed = 0ul;
    uint_t cre_passed = 0ul;
    for (uint_t iins = 0ul; iins < m_nop_insert; ++iins) {
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
