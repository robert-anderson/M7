//
// Created by rja on 28/07/22.
//

#ifndef M7_BOSHUBBARDUNIFORM_H
#define M7_BOSHUBBARDUNIFORM_H

#include "BosExcitGen.h"

namespace exgen {
    struct BosHubbardUniform : BosExcitGen {
        BosHubbardUniform(const BosHam& h, PRNG& prng) :
                BosExcitGen(h, prng, {opsig::c_0011}, "uniform Bose-Hubbard hopping") {}

        bool draw_bos(OpSig exsig, const field::BosOnv& src, prob_t& prob, conn::BosOnv& conn) override {
            DEBUG_ONLY(exsig);
            DEBUG_ASSERT_EQ(exsig, opsig::c_0011, "this excitation generator is only suitable for exsig 0011");
            const auto &occs = src.m_decoded.m_occ_modes.get();
            const auto nconn_lcm = m_h.m_basis.m_lattice->m_lcm_le_nadj_max;
            prob = 1/double(occs.size());
            auto rand = m_prng.draw_uint(occs.size() * nconn_lcm);
            const auto imode_ann = occs[rand / nconn_lcm];
            conn.clear();
            conn.m_ann.add(imode_ann, 1ul);
            const auto nadj = m_h.m_basis.m_lattice->m_sparse_adj.nentry(imode_ann);
            DEBUG_ASSERT_FALSE(nconn_lcm%nadj,
                               "LCM of connections should be exactly divisible by number of adjacent modes");
            // get iterator to beginning of adjacent sites
            auto it_cre = m_h.m_basis.m_lattice->m_sparse_adj.cbegin(imode_ann);
            // move the iterator by the drawn offset
            std::advance(it_cre, rand % nadj);
            // dereference the site (mode) index from the final iterator
            conn.m_cre.add(it_cre->m_i, 1ul);
            // scale probability to account for the number of adjacent sites we could have picked
            prob /= nadj;
            return true;
        }

        prob_t prob_bos(const field::BosOnv& src, const conn::BosOnv& conn) const override {
            const auto &occs = src.m_decoded.m_occ_modes.get();
            return 1/double(occs.size()*m_h.m_basis.m_lattice->m_sparse_adj.nentry(conn.m_ann[0].m_imode));
        }

        uint_t approx_nconn(OpSig /*exsig*/, sys::Particles particles) const override {
            return particles.m_bos;
        }
    };
}

#endif //M7_BOSHUBBARDUNIFORM_H
