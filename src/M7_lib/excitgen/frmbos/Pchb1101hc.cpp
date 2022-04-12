//
// Created by rja on 05/04/2022.
//

#include "Pchb1101hc.h"

Pchb1101hc::Pchb1101hc(const FrmBosHam &h, PRNG &prng) :
        FrmBosExcitGen(h, prng, {exsig_utils::ex_1101, exsig_utils::ex_1110}, "precomputed heatbath"),
        m_pick_n_given_pq(utils::pow<2>(h.m_bd.m_frm.m_nspinorb), h.m_bd.m_bos.m_nmode) {
    const auto nmode = m_h.m_bd.m_bos.m_nmode;
    const auto nspinorb = m_h.m_bd.m_frm.m_nspinorb;
    std::vector<defs::prob_t> weights(nmode, 0.0);
    size_t pq = 0ul;
    log::info("Initializing pre-computed samplers for fermion-boson kinetic term...");
    if (mpi::on_node_i_am_root()) {
        for (size_t p = 0ul; p < nspinorb; ++p) {
            for (size_t q = 0ul; q < nspinorb; ++q) {
                weights.assign(nmode, 0.0);
                if (p!=q) {
                    for (size_t n = 0ul; n < nmode; ++n) {
                        auto element = m_h.get_coeff_1101(n, p, q);
                        weights[n] = std::abs(element);
                    }
                }
                m_pick_n_given_pq.update(pq, weights);
                ++pq;
            }
        }
        DEBUG_ASSERT_EQ(pq, nspinorb*nspinorb, "not all single fermion excitations handled");
    }
    mpi::barrier();
}

bool Pchb1101hc::draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                             defs::prob_t &prob, conn::FrmBosOnv &conn) {
    /*
     * draw random occupied and vacant fermion indices
     */

    if (!UniformSingles::draw_spin_conserve_fn(m_prng, src.m_frm, prob, conn.m_frm)) return false;
    DEBUG_ASSERT_NE(prob, 0.0, "non-null generated with non-zero probability");
    /*
     * the fermion-boson coupling V_npq is defined with precedence to the "bosonic excitation" case (exsig 1110)
     */
    const bool cre = exsig_utils::decode_nbos_cre(exsig);
    const auto p = cre ? conn.m_frm.m_cre[0]: conn.m_frm.m_ann[0];
    const auto q = cre ? conn.m_frm.m_ann[0]: conn.m_frm.m_cre[0];
    const auto pq = p*m_h.m_bd.m_frm.m_nspinorb+q;
    auto n = m_pick_n_given_pq.draw(pq, m_prng);
    conn.m_bos.clear();
    if (cre) {
        if (src.m_bos[n] == m_h.m_bd.m_bos.m_nboson_max) return false;
        conn.m_bos.m_cre.set(n);
    }
    else {
        if (src.m_bos[n] == 0ul) return false;
        conn.m_bos.m_ann.set(n);
    };
    prob *= std::abs(m_h.get_coeff_1101(n, p, q)) / (m_pick_n_given_pq.norm(pq));
    return true;
}
