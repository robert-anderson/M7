//
// Created by Robert J. Anderson on 04/04/2022.
//

#include "Pchb2200.h"

Pchb2200::Pchb2200(const FrmHam &h, PRNG &prng):
        FrmExcitGen(h, prng, {exsig::ex_double}, "precomputed heat-bath fermion doubles"),
        m_nspinorb_pair(m_h.m_basis.m_nspinorb_pair),
        m_pick_ab_given_ij(m_nspinorb_pair, m_nspinorb_pair) {
    std::vector<defs::prob_t> weights(m_nspinorb_pair, 0.0);
    size_t ij = 0ul;
    log::info("Initializing pre-computed heat bath sampling weights for doubles...");
    const auto nspinorb = m_h.m_basis.m_nspinorb;
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

bool Pchb2200::draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                          defs::ham_t &helem, conn::FrmOnv &conn) {
    DEBUG_ASSERT_EQ(exsig, exsig::ex_double, "this excitation generator is only suitable for exsig 2200");
    size_t i, j, a, b;
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto npair_elec = integer::nspair(occs.size());
    size_t ij = m_prng.draw_uint(npair_elec);
    integer::inv_strigmap(j, i, ij);
    // i and j are positions in the occ list, convert to orb inds_t:
    i = occs[i];
    j = occs[j];

    DEBUG_ASSERT_LT(i, j, "drawn indices should be strictly ordered");

    // re-pack the selected indices into the flat index of the aliser row
    ij = integer::strigmap(j, i); // i and j are orbital indices
    if (fptol::numeric_zero(m_pick_ab_given_ij.norm(ij))) {
        // can't have a valid excitation if the row norm is zero
        return false;
    }

    size_t ab = m_pick_ab_given_ij.draw(ij, m_prng);
    integer::inv_strigmap(b, a, ab); // a and b are spin orbital indices
    //ASSERT(i!=a && i!=b && j!=a && j!=b)

    // reject the excitation if the creation indices are already occupied
    if (src.get(a)) return false;
    if (src.get(b)) return false;

    conn.m_ann.set(i, j);
    conn.m_cre.set(a, b);
    helem = m_h.get_element_2200(src, conn);
    prob = std::abs(helem) / (m_pick_ab_given_ij.norm(ij) * npair_elec);
    DEBUG_ASSERT_LE(prob, 1.0, "excitation drawn with invalid probability");
    return !fptol::numeric_zero(prob);
}

bool Pchb2200::draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
    /*
     * need the helement to compute the probability so if it isn't actually needed, just dispose of it
     * in contrast to the generic case where it is not assumed that the helement must be computed to get the prob,
     * so the delegation between the two virtual methods is reversed with respect to the generic case.
     */
    defs::ham_t helem;
    return draw_h_frm(exsig, src, prob, helem, conn);
}

defs::prob_t Pchb2200::prob_h_frm(const field::FrmOnv &src, const conn::FrmOnv &conn, defs::ham_t helem) const {
    const auto &occs = src.m_decoded.m_simple_occs.get();
    const auto npair_elec = integer::nspair(occs.size());
    auto ij = integer::strigmap(conn.m_ann[1], conn.m_ann[0]);
    return std::abs(helem) / (m_pick_ab_given_ij.norm(ij)*npair_elec);
}

defs::prob_t Pchb2200::prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const {
    return prob_h_frm(src, conn, m_h.get_element_2200(src, conn));
}

size_t Pchb2200::approx_nconn(size_t exsig, sys::Particles particles) const {
    return particles.m_frm.m_npair * m_nspinorb_pair;
}
