//
// Created by rja on 26/08/2021.
//

#include "LadderHoppingUniform.h"

LadderHoppingUniform::LadderHoppingUniform(const Hamiltonian &h, PRNG &prng) :
        LadderExcitGen(h, prng, {exsig_utils::ex_1110, exsig_utils::ex_1101}), m_singles(h, prng) {}

bool LadderHoppingUniform::draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                conn::FrmBosOnv &conn) {
    /*
     * draw random occupied and vacant fermion indices
     */
    if (!m_singles.draw(exsig_utils::ex_single, src, orbs, prob, conn)) return false;
    DEBUG_ASSERT_NE(prob, 0.0, "non-null generated with non-zero probability");
    /*
     * draw a random boson mode
     */
    const bool cre = exsig_utils::decode_nbos_cre(exsig);
    auto n = m_prng.draw_uint(m_h.nsite());
    conn.m_bos.clear();
    if (cre) {
        if (src.m_bos[n] == m_h.m_nboson_max) return false;
        conn.m_bos.m_cre.set({n, 1});
    }
    else {
        if (src.m_bos[n] == 0ul) return false;
        conn.m_bos.m_ann.set({n, 1});
    };
    prob /= m_h.nsite();
    return true;
}
