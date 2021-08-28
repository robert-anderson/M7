//
// Created by rja on 26/08/2021.
//

#include "LadderHoppingPc.h"

LadderHoppingPc::LadderHoppingPc(const Hamiltonian &h, PRNG &prng) :
        LadderHoppingUniform(h, prng),
        m_pick_n_given_pq(h.nsite()*h.nsite(), h.nsite()) {
    const auto nmode = h.nsite();
    std::vector<defs::prob_t> weights(nmode, 0.0);
    size_t pq = 0ul;
    log::info("Initializing pre-computed samplers for fermion-boson kinetic term...");
    if (mpi::on_node_i_am_root()) {
        for (size_t p = 0ul; p < nmode; ++p) {
            for (size_t q = 0ul; q < nmode; ++q) {
                weights.assign(nmode, 0.0);
                if (p!=q) {
                    for (size_t n = 0ul; n < nmode; ++n) {
                        auto element = m_h.m_ladder.v(n, p, q);
                        weights[n] = std::abs(element);
                    }
                }
                m_pick_n_given_pq.update(pq, weights);
                ++pq;
            }
        }
        DEBUG_ASSERT_EQ(pq, nmode*nmode, "not all single fermion excitations handled");
    }
    mpi::barrier();
}

bool LadderHoppingPc::draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                         conn::FrmBosOnv &conn) {
    /*
     * draw random occupied and vacant fermion indices
     */
    if (!m_singles.draw(exsig_utils::ex_single, src, orbs, prob, conn)) return false;
    DEBUG_ASSERT_NE(prob, 0.0, "non-null generated with non-zero probability");
    /*
     * the fermion-boson coupling V_npq is defined with precedence to the "bosonic excitation" case (exsig 1110)
     */
    const bool cre = exsig_utils::decode_nbos_cre(exsig);
    const auto p = (cre ? conn.m_frm.m_cre[0]: conn.m_frm.m_ann[0])%src.nsite();
    const auto q = (cre ? conn.m_frm.m_ann[0]: conn.m_frm.m_cre[0])%src.nsite();
    const auto pq = p*m_h.nsite()+q;
    auto n = m_pick_n_given_pq.draw(pq, m_prng);
    conn.m_bos.clear();
    if (cre) {
        if (src.m_bos[n] == m_h.m_nboson_max) return false;
        conn.m_bos.m_cre.set({n, 1});
    }
    else {
        if (src.m_bos[n] == 0ul) return false;
        conn.m_bos.m_ann.set({n, 1});
    };
    prob *= std::abs(m_h.m_ladder.v(n, p, q)) / (m_pick_n_given_pq.norm(pq));
    return true;
}