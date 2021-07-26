//
// Created by rja on 02/05/2021.
//

#include "HubbardSingles.h"

bool HubbardSingles::draw(const fields::FrmOnv &src_onv, fields::FrmOnv &dst_onv, const OccupiedOrbitals &occs,
                          const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) {
    auto rand = m_prng.draw_uint(2 * occs.size());
    bool choose_left = rand / occs.size();
    auto occ = occs[choose_left ? rand - occs.size() : rand];
    auto left = m_pbc ? conn_utils::left_pbc(occ, src_onv.m_nsite) : conn_utils::left(occ, src_onv.m_nsite);
    auto right = m_pbc ? conn_utils::right_pbc(occ, src_onv.m_nsite) : conn_utils::right(occ, src_onv.m_nsite);

    prob = 0.5 / occs.size();
    size_t vac = ~0ul;
    vac = choose_left ? left : right;
    if (vac == ~0ul) return false;
    if (src_onv.get(vac)) return false;

    conn.clear();
    conn.add(occ, vac);
    conn.apply(src_onv, dst_onv);
    helem = m_h->m_frm.get_element_1(src_onv, conn);
    return !consts::float_nearly_zero(helem, 1e-12);
}
