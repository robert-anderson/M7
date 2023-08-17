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
void do_annihilation(const uintv_t &ispinorbs, const field::FrmOnv& mbf) {
    auto& ref = const_cast<field::FrmOnv&>(mbf);
    for (auto& i: ispinorbs) ref.clr(i);
}

/**
 * helper to undo the above so as to leave the MBF unchanged
 */
void undo_annihilation(const uintv_t &ispinorbs, const field::FrmOnv& mbf) {
    auto& ref = const_cast<field::FrmOnv&>(mbf);
    for (auto& i: ispinorbs) ref.set(i);
}

void NotfMaeFiller::refresh_ann_map(const uintv_t &ann_ispinorbs, const v_t<uintp_t> &ann_siv, field::RdmInds& rdm_inds) {
    /*
     * don't refresh unless the outer indices have changed. in the first usage all spinorbs will be zero so if the
     * largest value is 0, it may be assumed that this is the first refresh
     */
    if (rdm_inds.m_frm.m_ann == ann_ispinorbs && ann_ispinorbs.back() != 0ul) return;
    m_ri_map.clear();
    auto fn = [&](uint_t imbf){
        m_hist.m_row.jump(imbf);
        auto has_phase = half_excit_phase(ann_ispinorbs, imbf);
        do_annihilation(ann_ispinorbs, m_hist.m_row.m_mbf);
        auto& ri_row = m_ri_map.insert(m_hist.m_row.m_mbf);
        undo_annihilation(ann_ispinorbs, m_hist.m_row.m_mbf);
        ri_row.m_weight[0] = m_hist.m_row.m_weight[0] * (has_phase ? -1.0 : 1.0);
    };
    bitset_isect::foreach_in_siv(ann_siv, fn);
    m_ri_map.remap_if_due();
    rdm_inds.m_frm.m_ann = ann_ispinorbs;
}

void NotfMaeFiller::probe_ann_map(const uintv_t &cre_ispinorbs, const v_t<uintp_t> &cre_siv, field::RdmInds& rdm_inds) {
    rdm_inds.m_frm.m_cre = cre_ispinorbs;
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

void NotfMaeFiller::resolve_identity(const uintv_t &ann_ispinorbs, const v_t<uintp_t> &ann_siv, const uintv_t &cre_ispinorbs,
                                     const v_t<uintp_t> &cre_siv, RdmInds &rdm_inds) {
    REQUIRE_TRUE_ALL(m_rdms, "RDMs object must be non-null");
    refresh_ann_map(ann_ispinorbs, ann_siv, rdm_inds);
    probe_ann_map(cre_ispinorbs, cre_siv, rdm_inds);
}

void NotfMaeFiller::fill() {
    auto fn = [&](const uintv_t& ao, const v_t<uintp_t>& ais, const uintv_t& co, const v_t<uintp_t>& cis, field::RdmInds& rdm_inds) {
        resolve_identity(ao, ais, co, cis, rdm_inds);
    };
    fill_foreach_set_pair(fn, m_rdms->all_ranksigs());
    // set the entire norm on the root rank
    if (mpi::i_am_root()) m_rdms->m_total_norm.m_local = get_norm();
}

wf_comp_t NotfMaeFiller::get_norm() const {
    auto& row = m_hist.m_row;
    wf_comp_t norm = 0.0;
    for (row.restart(); row; ++row) norm += math::pow<2>(std::abs(row.m_weight[0]));
    return norm;
}

NotfMaeFiller::NotfMaeFiller(const Table<MbfWeightRow> &hist, Rdms *rdms) :
        m_hist(hist), m_rdms(rdms),
        m_ind_displ(mpi::evenly_shared_displ(m_hist.nrow_in_use())),
        m_ind_count(mpi::evenly_shared_count(m_hist.nrow_in_use())),
        m_occ_bitsets(make_occ_bitsets()),
        m_partial_occ_bitsets(make_occ_bitsets(m_ind_displ, m_ind_count)),
        m_occ_sivs(bitset_isect::bitset_to_siv_many(m_occ_bitsets)),
        m_partial_occ_sivs(bitset_isect::bitset_to_siv_many(m_partial_occ_bitsets)),
        m_work_conn(m_hist.m_row.m_mbf.m_basis), m_ri_map(m_hist.m_row) {
}

void NotfMaeFiller::fill(const Table<MbfWeightRow> &hist, Rdms *rdms) {
    NotfMaeFiller filler(hist, rdms);
    filler.fill();
}

bool NotfMaeFiller::half_excit_phase(const uintv_t &ispinorbs, const Mbf &mbf) {
    m_work_conn.m_ann.clear();
    for (auto& i: ispinorbs) m_work_conn.m_ann.add(i);
    return m_work_conn.phase(mbf);
}

bool NotfMaeFiller::half_excit_phase(const uintv_t &ispinorbs, uint_t ihist_mbf) {
    m_hist.m_row.jump(ihist_mbf);
    const auto& mbf = m_hist.m_row.m_mbf;
    return half_excit_phase(ispinorbs, mbf);
}
