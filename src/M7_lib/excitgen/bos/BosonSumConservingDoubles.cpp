//
// Created by Robert J. Anderson on 25/11/2021.
//

#include "BosonSumConservingDoubles.h"


BosonSumConservingDoubles::BosonSumConservingDoubles(const BosHam &h, PRNG &prng) :
        BosExcitGen(h, prng, {exsig::ex_0022}, "boson mode index sum conserving") {}

void BosonSumConservingDoubles::set_a_range(uint_t i, uint_t j, uint_t &min, uint_t &max, uint_t &nexclude) const {
    const auto nmode = m_h.m_basis.m_nmode;
    min = 0ul;
    // e.g. nmode = 5, i = 3, j = 3. a can't be 0 or 1, since b would need to be 6 or 5 respectively to conserve sum
    if (i + j >= nmode) min = 1 + i + j - nmode;
    // e.g. nmode = 5, i = 1, j = 2. a can't be 4 since b would need to be -1 to conserve sum
    max = std::min(i + j + 1, nmode);
    nexclude = 1 + (i != j);
}

uint_t BosonSumConservingDoubles::na(uint_t i, uint_t j) const {
    uint_t min, max, nexclude;
    set_a_range(i, j, min, max, nexclude);
    return (max - min) - nexclude;
}

bool BosonSumConservingDoubles::draw_bos(uint_t, const field::BosOnv &src, prob_t &prob, conn::BosOnv &conn) {
    const auto &op_inds = src.m_decoded.m_expanded.get();
    const auto nboson_pair = integer::nspair(op_inds.size());
    uint_t ij = m_prng.draw_uint(nboson_pair);
    uint_t i, j;
    integer::inv_strigmap(j, i, ij);
    // i and j are positions in the occ list, convert to orb uintv_t:
    DEBUG_ASSERT_LT(i, j, "picked i, j should be in ascending order");
    DEBUG_ASSERT_LT(i, op_inds.size(), "i index OOB");
    DEBUG_ASSERT_LT(j, op_inds.size(), "j index OOB");
    i = op_inds[i];
    j = op_inds[j];
    DEBUG_ASSERT_LT(i, src.m_nelement, "i index OOB");
    DEBUG_ASSERT_LT(j, src.m_nelement, "j index OOB");

    uint_t min, max, nexclude;
    set_a_range(i, j, min, max, nexclude);

    if (min + nexclude == max) return false;
    uint_t a = m_prng.draw_uint(min, max - nexclude);
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
    DEBUG_ASSERT_LT(a, src.m_nelement, "boson index OOB");
    DEBUG_ASSERT_LT(b, src.m_nelement, "boson index OOB");
    DEBUG_ASSERT_EQ(i + j, a + b, "sum should be conserved");
    conn.m_ann.set(i, j);
    conn.m_cre.set(a, b);
    return true;
}

prob_t BosonSumConservingDoubles::prob_bos(const field::BosOnv &src, const conn::BosOnv &conn) const {
    const auto nboson_pair = integer::nspair(src.m_decoded.m_expanded.get().size());
    const auto i = conn.m_ann[0].m_imode;
    const auto j = conn.m_ann[0].m_imode==1 ? conn.m_ann[1].m_imode : i;
    prob_t prob;
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
