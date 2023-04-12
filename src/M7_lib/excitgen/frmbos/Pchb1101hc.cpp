//
// Created by Robert J. Anderson on 05/04/2022.
//

#include "Pchb1101hc.h"
#include "M7_lib/util/Math.h"

exgen::Pchb1101hc::Pchb1101hc(const FrmBosHam& h, PRNG& prng) :
        FrmBosExcitGen(h, prng, {opsig::c_1101, opsig::c_1110}, "precomputed heatbath"),
        m_pick_n_given_pq(math::pow<2>(h.m_basis.m_frm.m_nspinorb), h.m_basis.m_bos.m_nmode) {
    const auto nmode = m_h.m_basis.m_bos.m_nmode;
    const auto nspinorb = m_h.m_basis.m_frm.m_nspinorb;
    v_t<prob_t> weights(nmode, 0.0);
    uint_t pq = 0ul;
    logging::info("Initializing pre-computed samplers for fermion-boson kinetic term...");
    if (mpi::on_node_i_am_root()) {
        for (uint_t p = 0ul; p < nspinorb; ++p) {
            for (uint_t q = 0ul; q < nspinorb; ++q) {
                weights.assign(nmode, 0.0);
                if (p!=q) {
                    for (uint_t n = 0ul; n < nmode; ++n) {
                        auto element = m_h.get_coeff_1110(n, p, q);
                        weights[n] = std::abs(element);
                    }
                }
                m_pick_n_given_pq.update_(pq, weights);
                ++pq;
            }
        }
        DEBUG_ASSERT_EQ(pq, nspinorb*nspinorb, "not all single fermion excitations handled");
    }
    mpi::barrier();
}

bool exgen::Pchb1101hc::draw_frmbos(OpSig exsig, const field::FrmBosOnv& src,
                             prob_t& prob, conn::FrmBosOnv& conn) {
    /*
     * draw random occupied and vacant fermion indices
     */

    if (!exgen::UniformSingles::draw_spin_conserve_fn(m_prng, src.m_frm, prob, conn.m_frm)) return false;
    DEBUG_ASSERT_NE(prob, 0.0, "non-null generated with non-zero probability");
    /*
     * the fermion-boson coupling V_npq is defined with precedence to the "bosonic excitation" case (exsig 1110)
     */
    const bool cre = exsig.nbos_cre();
    const auto p = cre ? conn.m_frm.m_cre[0]: conn.m_frm.m_ann[0];
    const auto q = cre ? conn.m_frm.m_ann[0]: conn.m_frm.m_cre[0];
    const auto pq = p*m_h.m_basis.m_frm.m_nspinorb+q;
    auto n = m_pick_n_given_pq.draw(pq, m_prng);
    conn.m_bos.clear();
    if (cre) {
        if (src.m_bos[n] == m_h.m_basis.m_bos.m_occ_cutoff) return false;
        conn.m_bos.m_cre.set(n);
    }
    else {
        if (src.m_bos[n] == 0ul) return false;
        conn.m_bos.m_ann.set(n);
    };
    prob *= std::abs(m_h.get_coeff_1110(n, p, q)) / (m_pick_n_given_pq.norm(pq));
    return true;
}

prob_t exgen::Pchb1101hc::prob_h_frmbos(const field::FrmBosOnv& src,
                                       const conn::FrmBosOnv& conn, ham_t helem) const {
    auto prob = exgen::UniformSingles::prob_spin_conserve_fn(src.m_frm, conn.m_frm);
    const bool cre = conn.exsig().nbos_cre();
    const auto p = cre ? conn.m_frm.m_cre[0]: conn.m_frm.m_ann[0];
    const auto q = cre ? conn.m_frm.m_ann[0]: conn.m_frm.m_cre[0];
    const auto pq = p*m_h.m_basis.m_frm.m_nspinorb+q;
    return prob * std::abs(helem) / (m_pick_n_given_pq.norm(pq));
}

prob_t exgen::Pchb1101hc::prob_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn) const {
    return prob_h_frmbos(src, conn, m_h.get_element(src, conn));
}

uint_t exgen::Pchb1101hc::approx_nconn(OpSig, sys::Particles) const {
    return m_h.m_basis.m_frm.m_nsite*m_h.m_basis.m_bos.m_nmode;
}