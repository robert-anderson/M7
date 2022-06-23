//
// Created by Robert J. Anderson on 25/11/2021.
//

#ifndef M7_BOSONSUMCONSERVINGDOUBLES_H
#define M7_BOSONSUMCONSERVINGDOUBLES_H

#include "BosExcitGen.h"

struct BosonSumConservingDoubles : public BosExcitGen {
    BosonSumConservingDoubles(const BosHam &h, PRNG &prng):
        BosExcitGen(h, prng, {exsig::ex_0022}, "boson mode index sum conserving") {}

private:
    /**
     * @param min
     *  minimum allowed value for the draw of a
     * @param max
     *  maximum allowed value for the draw of a
     * @param nexclude
     *  index a should not be assigned the same value as i or j, step over these
     */
    void set_a_range(size_t i, size_t j, size_t& min, size_t &max, size_t& nexclude) const{
        const auto nmode = m_h.m_basis.m_nmode;
        min = 0ul;
        // e.g. nmode = 5, i = 3, j = 3. a can't be 0 or 1, since b would need to be 6 or 5 respectively to conserve sum
        if (i + j >= nmode) min = 1 + i + j - nmode;
        // e.g. nmode = 5, i = 1, j = 2. a can't be 4 since b would need to be -1 to conserve sum
        max = std::min(i + j + 1, nmode);
        nexclude = 1 + (i != j);
    }

    /**
     * @return
     *  number of possible values for index a
     */
    size_t na(size_t i, size_t j) const {
        size_t min, max, nexclude;
        set_a_range(i, j, min, max, nexclude);
        return (max - min) - nexclude;
    }

public:

    bool draw_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) override {
        const auto &op_inds = src.m_decoded.m_expanded.get();
        const auto nmode = src.m_nelement;
        const auto nboson_pair = integer::nspair(op_inds.size());
        size_t ij = m_prng.draw_uint(nboson_pair);
        size_t i, j;
        integer::inv_strigmap(j, i, ij);
        // i and j are positions in the occ list, convert to orb inds_t:
        DEBUG_ASSERT_LT(i, j, "picked i, j should be in ascending order");
        DEBUG_ASSERT_LT(i, op_inds.size(), "i index OOB");
        DEBUG_ASSERT_LT(j, op_inds.size(), "j index OOB");
        i = op_inds[i];
        j = op_inds[j];
        DEBUG_ASSERT_LT(i, src.m_nelement, "i index OOB");
        DEBUG_ASSERT_LT(j, src.m_nelement, "j index OOB");

        size_t min, max, nexclude;
        set_a_range(i, j, min, max, nexclude);

        if (min + nexclude == max) return false;
        size_t a = m_prng.draw_uint(min, max - nexclude);
        // skip mode indices to prevent selection of i or j as creation operators
        if (a >= i)++a;
        if (i != j && a >= j)++a;
        DEBUG_ASSERT_LT(a, max, "a is OOB");

        if (i == j) {
            prob = src[i] * (src[i] - 1) / double(2 * nboson_pair);
        } else {
            prob = src[i] * src[j] / double(nboson_pair);
        }
        auto b = (i + j) - a;
        if (a > b)std::swap(a, b);
        if (a != b) prob *= 2;
        prob /= (max - min) - nexclude;
        DEBUG_ASSERT_TRUE(i != a && j != a && i != b && j != b,
                          "should never draw same index for a pair of creation and annihilation ops");
        DEBUG_ASSERT_LT(a, nmode, "boson index OOB");
        DEBUG_ASSERT_LT(b, nmode, "boson index OOB");
        DEBUG_ASSERT_EQ(i + j, a + b, "sum should be conserved");
        conn.m_ann.set(i, j);
        conn.m_cre.set(a, b);
    }

    defs::prob_t prob_bos(const field::BosOnv &src, const conn::BosOnv &conn) const override {
        const auto nboson_pair = integer::nspair(src.m_decoded.m_expanded.get().size());
        const auto i = conn.m_ann[0].m_imode;
        const auto j = conn.m_ann[0].m_imode==1 ? conn.m_ann[1].m_imode : i;
        defs::prob_t prob;
        if (i == j) {
            prob = src[i] * (src[i] - 1) / double(2 * nboson_pair);
        } else {
            prob = src[i] * src[j] / double(nboson_pair);
        }
        const auto a = conn.m_cre[0].m_imode;
        auto b = (i + j) - a;
        if (a != b) prob *= 2;
        return prob / na(i, j);
    }

};


#endif //M7_BOSONSUMCONSERVINGDOUBLES_H
