//
// Created by rja on 09/05/2020.
//

#include "HeatBathDoubles.h"

HeatBathDoubles::HeatBathDoubles(const FermionHamiltonian *h, PRNG &prng) :
        FermionExcitationGenerator(h, prng, 2), m_pick_ab_given_ij(m_norb_pair, m_norb_pair) {
    std::vector<defs::prob_t> weights(m_norb_pair, 0.0);
    size_t ij = 0ul;
    size_t ab = 0ul;
    std::cout << "Initializing pre-computed heat bath sampling weights for doubles..." << std::endl;
    if (mpi::on_node_i_am_root()) {
        for (size_t i = 0ul; i < m_nintind; ++i) {
            for (size_t j = 0ul; j < i; ++j) {
                weights.assign(m_norb_pair, 0.0);
                ab = 0ul;
                for (size_t a = 0ul; a < m_nintind; ++a) {
                    for (size_t b = 0ul; b < a; ++b) {
                        //if (a!=i && a!=j && b!=i && b!=j) { !TODO why does this restriction fail?
                        weights[ab] = std::abs(m_h->get_element_2(i, j, a, b));
                        //}
                        ++ab;
                    }
                }
                m_pick_ab_given_ij.update(ij, weights);
                ASSERT(!consts::float_is_zero(m_pick_ab_given_ij.norm(ij)))
                ++ij;
            }
        }
        ASSERT(ij == m_norb_pair)
        ASSERT(ab == m_norb_pair)
    }
    mpi::barrier();
#ifndef NDEBUG
    for (ij = 0ul; ij < m_norb_pair; ++ij) {
        ASSERT(m_pick_ab_given_ij.nprob() == m_norb_pair)
        ASSERT(!consts::float_is_zero(m_pick_ab_given_ij.norm(ij)))
    }
#endif
}

bool HeatBathDoubles::draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv,
                                 const OccupiedOrbitals &occ, const VacantOrbitals &vac, defs::prob_t &prob,
                                 defs::ham_t &helem, conn::AsFermionOnv &anticonn) {
    // just draw uniform ij TODO! int weighted ij
    // return false if invalid excitation generated, true otherwise
    size_t i, j, a, b;
    size_t ij = m_prng.draw_uint(m_nelec_pair);
    integer_utils::inv_strigmap(j, i, ij);
    // i and j are positions in the occ list, convert to orb inds:
    i = occ.m_inds[i];
    j = occ.m_inds[j];
    ASSERT(std::any_of(occ.m_inds.cbegin(), occ.m_inds.cbegin() + occ.m_nind,
                       [&i](const size_t &k) { return k == i; }));
    ASSERT(std::any_of(occ.m_inds.cbegin(), occ.m_inds.cbegin() + occ.m_nind,
                       [&j](const size_t &k) { return k == j; }));
    ASSERT(i < j);

    ij = integer_utils::strigmap(j, i); // i and j are orbital indices
    size_t ab = m_pick_ab_given_ij.draw(ij, m_prng);
    integer_utils::inv_strigmap(b, a, ab); // a and b are orbital indices
    //ASSERT(i!=a && i!=b && j!=a && j!=b)

    auto either_vac_in_array = [&a, &b](const size_t &k) { return k == a || k == b; };

    if (std::any_of(occ.m_inds.cbegin(), occ.m_inds.cbegin() + occ.m_nind, either_vac_in_array)) {
        return 0;
    }
    anticonn.zero();
    anticonn.add(i, j, a, b);
    anticonn.apply(src_fonv, dst_fonv);
    helem = m_h->get_element_2(anticonn);
    prob = std::abs(helem) / (m_pick_ab_given_ij.norm(ij) * m_nelec_pair);
    ASSERT(prob <= 1)
    if (consts::float_nearly_zero(prob, 1e-14)) {
        return false;
    }
    return true;
}