//
// Created by anderson on 09/08/2023.
//

#include "NotfMaeFiller.h"

v_t<uintv_t> NotfMaeFiller::make_occ_bitsets(uint_t displ, uint_t count) const {
    v_t<uintv_t> bitsets;
    const auto nspinorb = m_hist.m_row.m_mbf.m_basis.m_nspinorb;
    for (uint_t ispinorb=0ul; ispinorb < nspinorb; ++ispinorb) {
        uintv_t inds;
        auto row = m_hist.m_row;
        // restrict iteration to chunk of histogrammed MBFs to which this rank is assigned
        for (row.restart(displ); row.in_range(displ + count); ++row) {
            // if spinorb occ append MBF index to this ispinorb's set of MBFs which have it occupied
            if (row.m_mbf.get(ispinorb)) inds.push_back(row.index());
        }
        // convert vector of set positions to a multiword bitset and append to vector of such bitsets
        bitsets.emplace_back(bitset_isect::make_bitset(inds, m_hist.nrow_in_use()));
    }
    return bitsets;
}

v_t<uintv_t> NotfMaeFiller::make_occ_bitsets() const {
    return make_occ_bitsets(0, m_hist.nrow_in_use());
}

/**
 * helper to temporarily violate const correctness so as to modify the MBF in place
 */
void do_annihilation(const uintv_t &ann_ispinorbs, const field::FrmOnv& mbf) {
    auto& ref = const_cast<field::FrmOnv&>(mbf);
    for (auto& i: ann_ispinorbs) ref.clr(i);
}

/**
 * helper to undo the above so as to leave the MBF unchanged
 */
void undo_annihilation(const uintv_t &ann_ispinorbs, const field::FrmOnv& mbf) {
    auto& ref = const_cast<field::FrmOnv&>(mbf);
    for (auto& i: ann_ispinorbs) ref.set(i);
}

void NotfMaeFiller::refresh_ann_map(const uintv_t &ann_ispinorbs, const v_t<uintp_t> &ann_siv, field::RdmInds& rdm_inds) {
    // don't refresh unless the outer indices have changed
    if (rdm_inds.m_frm.m_ann == ann_ispinorbs) return;
    auto fn = [&](uint_t imbf){
        m_hist.m_row.jump(imbf);
        auto has_phase = half_excit_phase(ann_ispinorbs, imbf);
        do_annihilation(ann_ispinorbs, m_hist.m_row.m_mbf);
        auto& ri_row = m_ri_map.insert(m_hist.m_row.m_mbf);
        undo_annihilation(ann_ispinorbs, m_hist.m_row.m_mbf);
        ri_row.m_weight[0] = m_hist.m_row.m_weight[0] * (has_phase ? -1.0 : 1.0);
    };
    bitset_isect::foreach_in_siv(ann_siv, fn);
    rdm_inds.m_frm.m_ann = ann_ispinorbs;
}

void NotfMaeFiller::probe_ann_map(const uintv_t &cre_ispinorbs, const v_t<uintp_t> &cre_siv, field::RdmInds& rdm_inds) {
    auto fn = [&](uint_t imbf){
        m_hist.m_row.jump(imbf);
        auto has_phase = half_excit_phase(cre_ispinorbs, imbf);
        do_annihilation(cre_ispinorbs, m_hist.m_row.m_mbf);
        auto& ri_row = m_ri_map.lookup(m_hist.m_row.m_mbf);
        undo_annihilation(cre_ispinorbs, m_hist.m_row.m_mbf);
        if (!ri_row) return;
        const auto contrib = ri_row.m_weight[0] * m_hist.m_row.m_weight[0] * (has_phase ? -1.0 : 1.0);
        m_rdms->make_full_contrib(rdm_inds, contrib, false);
    };
    bitset_isect::foreach_in_siv(cre_siv, fn);
}