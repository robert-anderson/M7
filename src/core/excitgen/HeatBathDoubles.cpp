//
// Created by rja on 09/05/2020.
//

#include "HeatBathDoubles.h"

HeatBathDoubles::HeatBathDoubles(const Hamiltonian &h, PRNG &prng) :
        FrmExcitGen(h, prng, 2), m_pick_ab_given_ij(m_norb_pair, m_norb_pair) {
    std::vector<defs::prob_t> weights(m_norb_pair, 0.0);
    size_t ij = 0ul;
    log::info("Initializing pre-computed heat bath sampling weights for doubles...");
    if (mpi::on_node_i_am_root()) {
        for (size_t i = 0ul; i < m_nspinorb; ++i) {
            for (size_t j = 0ul; j < i; ++j) {
                weights.assign(m_norb_pair, 0.0);
                size_t ab = 0ul;
                for (size_t a = 0ul; a < m_nspinorb; ++a) {
                    for (size_t b = 0ul; b < a; ++b) {
                        //if (a!=i && a!=j && b!=i && b!=j) { !TODO why does this restriction fail?
                        auto element = m_h.m_frm.get_element_2200(i, j, a, b);
                        weights[ab] = std::abs(element);
                        //}
                        ++ab;
                    }
                }
                ASSERT(ab == m_norb_pair)
                m_pick_ab_given_ij.update(ij, weights);
                ++ij;
            }
        }
        ASSERT(ij == m_norb_pair)
    }
    mpi::barrier();
#ifndef NDEBUG
    for (ij = 0ul; ij < m_norb_pair; ++ij) {
        ASSERT(m_pick_ab_given_ij.nprob() == m_norb_pair)
    }
#endif
}

size_t HeatBathDoubles::approx_nconn() const {
    return m_nelec_pair*m_norb_pair;
}
