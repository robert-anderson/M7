//
// Created by rja on 02/05/2021.
//

#include "HubbardSingles.h"

bool HubbardSingles::_draw(const fields::FermionOnv &src_onv, fields::FermionOnv &dst_onv, const OccupiedOrbitals &occs,
                           const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
                           conn::Antisym<0> &anticonn) {
    auto rand = m_prng.draw_uint(2*occs.size());
    bool choose_left = rand/occs.size();
    auto occ = occs[choose_left ? rand-occs.size() : rand];
    auto left = conn_utils::left(occ, src_onv.m_nsite);
    auto right = conn_utils::right(occ, src_onv.m_nsite);

    prob = 0.5/occs.size();
    size_t vac = ~0ul;
    vac = choose_left ? left : right;
    if (vac==~0ul) return false;
    if (src_onv.get(vac)) return false;

    anticonn.zero();
    anticonn.add(occ, vac);
    anticonn.apply(src_onv, dst_onv);
    helem = m_h->get_element_1(anticonn);
    return !consts::float_nearly_zero(helem, 1e-12);
}
