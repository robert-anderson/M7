//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHDOUBLES_H
#define M7_HEATBATHDOUBLES_H

#include "ExcitGen.h"
#include "src/core/sample/Aliaser.h"
#include "src/core/field/Fields.h"

/**
 * precomputed sampler for fermion double excitations
 */
class HeatBathDoubles : public FrmExcitGen {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const Hamiltonian &h, PRNG &prng);

    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override {
        // just draw uniform ij TODO! int weighted ij
        // return false if invalid excitation generated, true otherwise
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

    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_frm, orbs, prob, helem, conn.m_frm);
    }

    size_t approx_nconn() const override;

};

#endif //M7_HEATBATHDOUBLES_H
