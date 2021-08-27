//
// Created by rja on 09/05/2020.
//

#include "HeatBathDoubles.h"

HeatBathDoubles::HeatBathDoubles(const Hamiltonian &h, PRNG &prng) :
        FrmExcitGen(h, prng, exsig_utils::ex_double), m_pick_ab_given_ij(m_norb_pair, m_norb_pair) {
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

bool HeatBathDoubles::draw(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                           defs::ham_t &helem, conn::FrmOnv &conn) {
    size_t i, j, a, b;
    size_t ij = m_prng.draw_uint(m_nelec_pair);
    integer_utils::inv_strigmap(j, i, ij);
    const auto& occs = orbs.occ(src).m_flat;
    // i and j are positions in the occ list, convert to orb inds:
    i = occs[i];
    j = occs[j];
    ASSERT(std::any_of(occs.inds().cbegin(), occs.inds().cend(),
                       [&i](const size_t &k) { return k == i; }));
    ASSERT(std::any_of(occs.inds().cbegin(), occs.inds().cend(),
                       [&j](const size_t &k) { return k == j; }));
    ASSERT(i < j);

    ij = integer_utils::strigmap(j, i); // i and j are orbital indices
    if (consts::float_is_zero(m_pick_ab_given_ij.norm(ij))){
        // can't have a valid excitation if the row norm is zero
        return false;
    }

    size_t ab = m_pick_ab_given_ij.draw(ij, m_prng);
    integer_utils::inv_strigmap(b, a, ab); // a and b are orbital indices
    //ASSERT(i!=a && i!=b && j!=a && j!=b)

    auto either_vac_in_array = [&a, &b](const size_t &k) { return k == a || k == b; };

    if (std::any_of(occs.inds().cbegin(), occs.inds().end(), either_vac_in_array)) {
        return false;
    }
    conn.set(i, j, a, b);
    helem = m_h.m_frm.get_element_2200(src, conn);
    prob = std::abs(helem) / (m_pick_ab_given_ij.norm(ij) * m_nelec_pair);
    DEBUG_ASSERT_LE(prob, 1.0, "excitation drawn with invalid probability");
    if (consts::float_nearly_zero(prob, 1e-14)) {
        return false;
    }
    return true;
}

size_t HeatBathDoubles::approx_nconn() const {
    return m_nelec_pair*m_norb_pair;
}