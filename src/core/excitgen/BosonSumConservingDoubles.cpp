//
// Created by rja on 25/11/2021.
//

#include "BosonSumConservingDoubles.h"

BosonSumConservingDoubles::BosonSumConservingDoubles(const Hamiltonian &h, PRNG &prng) :
        BosExcitGen(h, prng, exsig_utils::ex_0022),
        m_nboson_pair(integer_utils::nspair(h.nboson())) {}

bool BosonSumConservingDoubles::draw(const size_t &exsig, const BosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                                     conn::BosOnv &conn) {
    const auto &op_inds = orbs.bos_op_inds(src);
    const auto &nmode = src.m_nelement;
    const auto &nboson = m_h.nboson();
    size_t ij = m_prng.draw_uint(m_nboson_pair);
    size_t i, j;
    integer_utils::inv_strigmap(j, i, ij);
    // i and j are positions in the occ list, convert to orb inds:
    i = op_inds[i];
    j = op_inds[j];
    size_t nexclude = 1 + (i != j);
    // a must be >= i+j-(nmode-1)
    // so that b (= i+j - a) is <= nmode-1
    // a must be <= min(i+j, nmode)
    auto min = 1 + i + j;
    if (min>nmode) min-=nmode;
    else min = 0;
    auto max = std::min(i + j, nmode) - nexclude;
    if (min == max) return false;
    size_t a = m_prng.draw_uint(min, max);
    if (a >= i) ++a;
    if (i != j && a >= j) ++a;
    DEBUG_ASSERT_LT(a, std::min(i + j, nmode), "b=i+j-a is not possible from this value of a");
    size_t b = (i + j) - a;
    if (a > b) std::swap(a, b);

    if (i == j) {
        prob = src[i] * (src[i]-1) / double(2*nboson * (nboson - 1));
    } else {
        prob = src[i] * src[j] / double(nboson * (nboson - 1));
    }
    if (a!=b) prob*=2;
    prob /= max - min;
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
