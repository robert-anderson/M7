//
// Created by rja on 04/04/2022.
//

#include "PchbDoubles2.h"

PchbDoubles2::PchbDoubles2(const FrmHam &h, PRNG &prng) :
        FrmExcitGen2(h, prng, {exsig_utils::ex_double}, "precomputed heat-bath fermion doubles"),
        m_pick_ab_given_ij(m_nspinorb_pair, m_nspinorb_pair) {
    std::vector<defs::prob_t> weights(m_nspinorb_pair, 0.0);
    size_t ij = 0ul;
    log::info("Initializing pre-computed heat bath sampling weights for doubles...");
    const auto nspinorb = 2 * m_h.m_nsite;
    if (mpi::on_node_i_am_root()) {
        for (size_t i = 0ul; i < nspinorb; ++i) {
            for (size_t j = 0ul; j < i; ++j) {
                weights.assign(m_nspinorb_pair, 0.0);
                size_t ab = 0ul;
                for (size_t a = 0ul; a < nspinorb; ++a) {
                    for (size_t b = 0ul; b < a; ++b) {
                        //if (a!=i && a!=j && b!=i && b!=j) { !TODO why does this restriction fail?
                        auto element = m_h.get_coeff_2200(i, j, a, b);
                        weights[ab] = std::abs(element);
                        //}
                        ++ab;
                    }
                }
                DEBUG_ASSERT_EQ(ab, m_nspinorb_pair, "loop did not include all vacant orbital pairs");
                m_pick_ab_given_ij.update(ij, weights);
                ++ij;
            }
        }
        DEBUG_ASSERT_EQ(ij, m_nspinorb_pair, "loop did not include all occupied orbital pairs");
    }
    mpi::barrier();
}

bool PchbDoubles2::draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                             defs::ham_t &helem, conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig_utils::ex_double, "this excitation generator is only suitable for exsig 2200");
    size_t i, j, a, b;
    size_t ij = m_prng.draw_uint(m_nelec_pair);
    integer_utils::inv_strigmap(j, i, ij);
    const auto &occs = src.m_decoded.m_simple_occs.get();
    DEBUG_ASSERT_EQ(occs.size(), m_h.m_nelec, "incorrect number of electrons")
    // i and j are positions in the occ list, convert to orb inds:
    i = occs[i];
    j = occs[j];

    DEBUG_ASSERT_LT(i, j, "drawn indices should be strictly ordered");

    // re-pack the selected indices into the flat index of the aliser row
    ij = integer_utils::strigmap(j, i); // i and j are orbital indices
    if (consts::nearly_zero(m_pick_ab_given_ij.norm(ij))) {
        // can't have a valid excitation if the row norm is zero
        return false;
    }

    size_t ab = m_pick_ab_given_ij.draw(ij, m_prng);
    integer_utils::inv_strigmap(b, a, ab); // a and b are spin orbital indices
    //ASSERT(i!=a && i!=b && j!=a && j!=b)

    // reject the excitation if the creation indices are already occupied
    if (src.get(a)) return false;
    if (src.get(b)) return false;

    conn.set(i, j, a, b);
    helem = m_h.get_element_2200(src, conn);
    prob = std::abs(helem) / (m_pick_ab_given_ij.norm(ij) * m_nelec_pair);
    DEBUG_ASSERT_LE(prob, 1.0, "excitation drawn with invalid probability");
    if (consts::nearly_zero(prob, 1e-14)) {
        return false;
    }
    return true;
}

bool PchbDoubles2::draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
    /*
     * need the helement to compute the probability so if it isn't actually needed, just dispose of it
     * in contrast to the generic case where it is not assumed that the helement must be computed to get the prob,
     * so the delegation between the two virtual methods is reversed with respect to the generic case.
     */
    defs::ham_t helem;
    return draw_h_frm(exsig, src, prob, helem, conn);
}

size_t PchbDoubles2::approx_nconn() const {
    return m_nelec_pair * m_nspinorb_pair;
}
