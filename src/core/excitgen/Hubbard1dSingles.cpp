//
// Created by rja on 02/05/2021.
//

#include "Hubbard1dSingles.h"

size_t Hubbard1dSingles::approx_nconn() const {
    return m_nelec;
}

bool Hubbard1dSingles::draw(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs,
                            defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) {
    const auto nelec = m_h.nelec();
    auto rand = m_prng.draw_uint(2 * nelec);
    bool choose_left = rand / nelec;
    const auto &occs = orbs.occ(src).m_flat;
    auto occ = occs[choose_left ? rand - occs.size() : rand];
    auto left = m_pbc ? model_utils::left_pbc(occ, src.m_nsite) : model_utils::left(occ, src.m_nsite);
    auto right = m_pbc ? model_utils::right_pbc(occ, src.m_nsite) : model_utils::right(occ, src.m_nsite);

    prob = 0.5 / occs.size();
    size_t vac = ~0ul;
    vac = choose_left ? left : right;
    if (vac == ~0ul) return false;
    if (src.get(vac)) return false;

    conn.set(occ, vac);
    helem = m_h.m_frm.get_element_1100(src, conn);
    return !consts::float_nearly_zero(helem, 1e-12);
}

bool Hubbard1dSingles::draw(const size_t &exsig, const FrmBosOnv &src_onv, CachedOrbs &orbs,
                            defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
    return draw(exsig, src_onv.m_frm, orbs, prob, helem, conn.m_frm);
}

std::string Hubbard1dSingles::description() const {
    return log::format("Hubbard 1D with {}BCs", m_pbc ? "P" : "O");
}
