//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/foreach/ForeachVirtual.h>
#include "src/core/caches/CachedOrbs.h"
#include "gtest/gtest.h"

TEST(DecodedDeterminant, SimpleOccAndVac){
    using namespace decoded_mbf::frm;
    buffered::FrmOnv mbf({50, 0});
    defs::inds setbits{0, 1, 4, 7, 32, 50, 51, 54, 60, 89, 99};
    mbf = setbits;
    defs::inds clrbits;
    auto iter = setbits.begin();
    for (size_t i=0ul; i < mbf.nbit(); ++i){
        if (iter!=setbits.end() && i==*iter) iter++;
        else clrbits.push_back(i);
    }

    auto& occs = mbf.m_decoded.m_simple_occs.get();
    ASSERT_TRUE(std::equal(occs.cbegin(), occs.cend(), setbits.cbegin()));

    auto& vacs = mbf.m_decoded.m_simple_vacs.get();
    ASSERT_TRUE(std::equal(vacs.cbegin(), vacs.cend(), clrbits.cbegin()));

    /*
     * for a small number of occupied orbs, run through all possible arrangements
     */
    const size_t noccorb = 3;
    auto occ_fn = [&mbf](size_t iiter, const defs::inds &value) {
        mbf.zero();
        mbf = value;
        mbf.m_decoded.clear();
        auto& occ_simple_inds = mbf.m_decoded.m_simple_occs.get();
        ASSERT_EQ(occ_simple_inds, value);
    };
    foreach_virtual::rtnd::lambda::Ordered<> occ_foreach(occ_fn, mbf.m_nspinorb, noccorb);
    occ_foreach.loop();

    /*
     * for a small number of vacant orbs, run through all possible arrangements
     */
    const size_t nvacorb = 3;
    auto vac_fn = [&mbf](size_t iiter, const defs::inds &value) {
        mbf.set();
        for (auto i: value) mbf.clr(i);
        mbf.m_decoded.clear();
        auto& vac_simple_inds = mbf.m_decoded.m_simple_vacs.get();
        ASSERT_EQ(vac_simple_inds, value);
    };
    foreach_virtual::rtnd::lambda::Ordered<> vac_foreach(vac_fn, mbf.m_nspinorb, nvacorb);
    vac_foreach.loop();
}

TEST(DecodedDeterminant, SpinSymOccAndVac){
    using namespace decoded_mbf::frm;
    // arbitrary, fictitious group
    AbelianGroup grp({"X", "Y", "Z"}, [](const size_t& iirrep, const size_t& jirrep){
        return (iirrep+jirrep)%3;
    });

    /*                                         0  1  2  3  4  5  6  7  8  9
     * occupied orbitals (alpha):              o  o  o  /  o  /  /  o  /  o
     * occupied orbitals (beta ):              /  /  o  /  o  o  o  o  o  /
     */
    AbelianGroupMap grp_map(grp, {0, 1, 1, 2, 0, 0, 1, 2, 1, 0});
    ASSERT_EQ(grp_map.m_nsite, 10);
    /*
     * irrep occupations (alpha):
     *  0("Xa"): 3, 1("Ya"): 2, 2("Za"): 1
     *
     * irrep vacancies (alpha):
     *  0("Xa"): 1, 1("Ya"): 2, 2("Za"): 1
     *
     *
     * irrep occupations (beta):
     *  3("Xb"): 2, 4("Yb"): 3, 5("Zb"): 1
     *
     * irrep vacancies (beta):
     *  3("Xb"): 2, 4("Yb"): 1, 5("Zb"): 1
     */

    buffered::FrmOnv mbf({10, 0, grp_map});
    defs::inds alpha_occ{0, 1, 2, 4, 7, 9};
    defs::inds beta_occ{2, 4, 5, 6, 7, 8};
    mbf = {alpha_occ, beta_occ};
    ASSERT_EQ(mbf.nsetbit(), alpha_occ.size() + beta_occ.size());

    defs::inds chk_inds;

    auto& occs = mbf.m_decoded.m_spin_sym_occs.get();
    auto& vacs = mbf.m_decoded.m_spin_sym_vacs.get();
    /*
     * alpha occupied
     */
    ASSERT_EQ(occs.size(0), 3);
    chk_inds = occs[{0, 0}];
    ASSERT_EQ(chk_inds, defs::inds({0, 4, 9}));
    ASSERT_EQ(occs.size(1), 2);
    chk_inds = occs[{0, 1}];
    ASSERT_EQ(chk_inds, defs::inds({1, 2}));
    ASSERT_EQ(occs.size(2), 1);
    chk_inds = occs[{0, 2}];
    ASSERT_EQ(chk_inds, defs::inds({7}));

    /*
     * alpha vacant
     */
    ASSERT_EQ(vacs.size(0), 1);
    chk_inds = vacs[{0, 0}];
    ASSERT_EQ(chk_inds, defs::inds({5}));
    ASSERT_EQ(vacs.size(1), 2);
    chk_inds = vacs[{0, 1}];
    ASSERT_EQ(chk_inds, defs::inds({6, 8}));
    ASSERT_EQ(vacs.size(2), 1);
    chk_inds = vacs[{0, 2}];
    ASSERT_EQ(chk_inds, defs::inds({3}));

    /*
     * beta occupied
     * beta_occ{2, 4, 5, 6, 7, 8} -> 12, 14, 15, 16, 17, 18
     */
    ASSERT_EQ(occs.size(3), 2);
    chk_inds = occs[{1, 0}];
    ASSERT_EQ(chk_inds, defs::inds({14, 15}));
    ASSERT_EQ(occs.size(4), 3);
    chk_inds = occs[{1, 1}];
    ASSERT_EQ(chk_inds, defs::inds({12, 16, 18}));
    ASSERT_EQ(occs.size(5), 1);
    chk_inds = occs[{1, 2}];
    ASSERT_EQ(chk_inds, defs::inds({17}));

    /*
     * beta vacant
     */
    ASSERT_EQ(vacs.size(3), 2);
    chk_inds = vacs[{1, 0}];
    ASSERT_EQ(chk_inds, defs::inds({10, 19}));
    ASSERT_EQ(vacs.size(4), 1);
    chk_inds = vacs[{1, 1}];
    ASSERT_EQ(chk_inds, defs::inds({11}));
    ASSERT_EQ(vacs.size(5), 1);
    chk_inds = vacs[{1, 2}];
    ASSERT_EQ(chk_inds, defs::inds({13}));
}