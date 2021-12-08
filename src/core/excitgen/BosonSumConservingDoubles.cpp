//
// Created by rja on 25/11/2021.
//

#include "BosonSumConservingDoubles.h"

BosonSumConservingDoubles::BosonSumConservingDoubles(const Hamiltonian &h, PRNG &prng) :
        BosExcitGen(h, prng, exsig_utils::ex_0022),
        m_nboson_pair(integer_utils::nspair(h.nboson())) {}

bool BosonSumConservingDoubles::draw_bos(const size_t &exsig, const BosOnv &src,
                                         CachedOrbs &orbs, defs::prob_t &prob, conn::BosOnv &conn) {
    const auto &op_inds = orbs.bos_op_inds(src);
    DEBUG_ASSERT_EQ(op_inds.size(), m_h.nboson(), "picked i, j should be in ascending order");
    const auto &nmode = src.m_nelement;
    size_t ij = m_prng.draw_uint(m_nboson_pair);
    size_t i, j;
    integer_utils::inv_strigmap(j, i, ij);
    // i and j are positions in the occ list, convert to orb inds:
    DEBUG_ASSERT_LT(i, j, "picked i, j should be in ascending order");
    DEBUG_ASSERT_LT(i, m_h.nboson(), "i index OOB");
    DEBUG_ASSERT_LT(j, m_h.nboson(), "j index OOB");
    i = op_inds[i];
    j = op_inds[j];
    DEBUG_ASSERT_LT(i, src.m_nelement, "i index OOB");
    DEBUG_ASSERT_LT(j, src.m_nelement, "j index OOB");

    auto min = 0ul;
    // e.g. nmode = 5, i = 3, j = 3. a can't be 0 or 1, since b would need to be 6 or 5 respectively to conserve sum
    if (i + j >= nmode) min = 1 + i + j - nmode;
    // e.g. nmode = 5, i = 1, j = 2. a can't be 4 since b would need to be -1 to conserve sum
    auto max = std::min(i + j + 1, nmode);

    size_t nexclude = 1 + (i != j);
    if (min + nexclude == max) return false;
    size_t a = m_prng.draw_uint(min, max - nexclude);
    // skip mode indices to prevent selection of i or j as creation operators
    if (a >= i)++a;
    if (i != j && a >= j)++a;
    DEBUG_ASSERT_LT(a, max, "a is OOB");

    if (i == j) {
        prob = src[i] * (src[i] - 1) / double(2 * m_nboson_pair);
    } else {
        prob = src[i] * src[j] / double(m_nboson_pair);
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
    return true;
}

std::string BosonSumConservingDoubles::description() const {
    return "Boson sum conserving";
}
