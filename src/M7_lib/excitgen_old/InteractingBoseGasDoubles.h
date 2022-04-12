//
// Created by rja on 07/03/2022.
//

#ifndef M7_INTERACTINGBOSEGASDOUBLES_H
#define M7_INTERACTINGBOSEGASDOUBLES_H

#include "M7_lib/hamiltonian/frmbos/InteractingBoseGasBosHam.h"

#include "ExcitGen.h"

class InteractingBoseGasDoubles : BosExcitGen {
    const size_t m_nbos_pair;
    /**
     * working arrays into which the mode indices are decoded
     */
    mutable std::vector<int> m_ikpoints_i_work, m_ikpoints_j_work;//, m_ikpoints_k_work, m_ikpoints_l_work;
    InteractingBoseGasDoubles(const Hamiltonian& h, PRNG& prng) :
        BosExcitGen(h, prng, exsig_utils::ex_0022),
        m_nbos_pair(integer_utils::combinatorial(h.nboson(), 2)){}

    const InteractingBoseGasBosHam* h_cast() const {
        return dynamic_cast<const InteractingBoseGasBosHam*>(m_h.m_bos.get());
    }

    bool draw_h_bos(const size_t &exsig, const BosOnv &src, CachedOrbs &orbs, defs::prob_t &prob, defs::ham_t &helem,
                    conn::BosOnv &conn) override {
        auto h = h_cast();
        DEBUG_ONLY(h);
        DEBUG_ASSERT_TRUE(h, "Boson hamiltonian couldn't be cast to InteractingBosGasBosHam");
        auto& bos_op_inds = orbs.bos_op_inds(src);
        DEBUG_ASSERT_EQ(bos_op_inds.size(), m_h.nboson(),
                        "boson ONV does not hold expected number of particles");
        auto ipair_occ = m_prng.draw_uint(m_nbos_pair);
        size_t imode, jmode;
        integer_utils::inv_strigmap(imode, jmode, ipair_occ);
        imode = bos_op_inds[imode];
        jmode = bos_op_inds[jmode];
//        h->get_ikpoints(imode, m_ikpoints_i_work);
//        h->get_ikpoints(jmode, m_ikpoints_j_work);
        // sum the momenta in each dim
        for (size_t i=0ul; i<m_ikpoints_i_work.size(); ++i) m_ikpoints_i_work[i]+=m_ikpoints_j_work[i];
        /*
         * ultimately, i+j = k+l must hold, so the picking of k from here, defines l
         * must pick
         * k = i+j-l
         * k >= i+j-l_min
         * k <= i+j-l_max
         */
        return true;
    }

    bool draw_bos(const size_t &exsig, const BosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                  conn::BosOnv &conn) override {
        defs::ham_t helem;
        return draw_h_bos(exsig, src, orbs, prob, helem, conn);
    }

};


#endif //M7_INTERACTINGBOSEGASDOUBLES_H
