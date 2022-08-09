//
// Created by Robert John Anderson on 2020-03-31.
//

#include <M7_lib/table/BufferedFields.h>
#include "gtest/gtest.h"

TEST(DecodedMbf, Simple){
    const uintv_t setbits{0, 1, 4, 7, 32, 50, 51, 54, 60, 89, 99};
    buffered::FrmOnv mbf(50);
    mbf = setbits;
    uintv_t clrbits;
    auto iter = setbits.begin();
    for (uint_t i=0ul; i < mbf.nbit(); ++i){
        if (iter!=setbits.end() && i==*iter) iter++;
        else clrbits.push_back(i);
    }

    auto& occs = mbf.m_decoded.m_simple_occs.get();
    ASSERT_TRUE(std::equal(occs.cbegin(), occs.cend(), setbits.cbegin()));

    auto& vacs = mbf.m_decoded.m_simple_vacs.get();
    ASSERT_TRUE(std::equal(vacs.cbegin(), vacs.cend(), clrbits.cbegin()));

    using namespace basic_foreach::rtnd;
    /*
     * for a small number of occupied orbs, run through all possible arrangements
     */
    const uint_t noccorb = 3;
    auto occ_fn = [&mbf](const inds_t& inds) {
        mbf.zero();
        mbf = inds;
        mbf.m_decoded.clear();
        auto& occ_simple_inds = mbf.m_decoded.m_simple_occs.get();
        ASSERT_EQ(occ_simple_inds, inds);
    };
    Ordered<> occ_foreach(mbf.m_basis.m_nspinorb, noccorb);
    occ_foreach.loop(occ_fn);

    /*
     * for a small number of vacant orbs, run through all possible arrangements
     */
    const uint_t nvacorb = 3;
    auto vac_fn = [&mbf](const uintv_t &inds) {
        mbf.set();
        for (auto i: inds) mbf.clr(i);
        mbf.m_decoded.clear();
        auto& vac_simple_inds = mbf.m_decoded.m_simple_vacs.get();
        ASSERT_EQ(vac_simple_inds, inds);
    };
    Ordered<> vac_foreach(mbf.m_basis.m_nspinorb, nvacorb);
    vac_foreach.loop(vac_fn);
}

TEST(DecodedMbf, SingleMultipleOccupation) {
    const uint_t nsite = 8;
    uintv_t alpha_inds = {0, 4, 5, 7};
    uintv_t beta_inds = {0, 1, 2, 5};
    buffered::FrmOnv mbf(nsite);
    mbf = {alpha_inds, beta_inds};
    uintv_t chk_inds;

    chk_inds = {0, 1, 2, 4, 5, 7};
    ASSERT_EQ(mbf.m_decoded.m_occ_sites.get(), chk_inds);

    chk_inds = {0, 5};
    ASSERT_EQ(mbf.m_decoded.m_doubly_occ_sites.get(), chk_inds);

    chk_inds = {0, 3, 5, 6};
    ASSERT_EQ(mbf.m_decoded.m_not_singly_occ_sites.get(), chk_inds);
}

TEST(DecodedMbf, Labelled){
    // arbitrary, fictitious group
    AbelianGroup grp({"W", "X", "Y", "Z"}, [](const uint_t& iirrep, const uint_t& jirrep){
        return (iirrep+jirrep)%4;
    });

    uintv_t alpha_occ{0, 1, 2, 4, 7, 9};
    uintv_t beta_occ{2, 4, 5, 6, 7, 8};
    /*                                         0  1  2  3  4  5  6  7  8  9
     * occupied orbitals (alpha):              o  o  o  /  o  /  /  o  /  o
     * occupied orbitals (beta ):              /  /  o  /  o  o  o  o  o  /
     */
    AbelianGroupMap grp_map(grp, {0, 1, 1, 3, 0, 0, 1, 3, 1, 0});
    // note that label Y is unused
    ASSERT_EQ(grp_map.m_nsite, 10);
    /*
     * irrep occupations (alpha):
     *  0("Wa"): 3, 1("Xa"): 2, 2("Ya"): 0, 3("Za"): 1
     *
     * irrep vacancies (alpha):
     *  0("Wa"): 1, 1("Xa"): 2, 2("Ya"): 0, 3("Za"): 1
     *
     *
     * irrep occupations (beta):
     *  4("Wb"): 2, 5("Xb"): 3, 6("Yb"): 0, 7("Zb"): 1
     *
     * irrep vacancies (beta):
     *  4("Wb"): 2, 5("Xb"): 1, 6("Yb"): 0, 7("Zb"): 1
     */

    buffered::FrmOnv mbf({10, grp_map});
    mbf = {alpha_occ, beta_occ};
    ASSERT_EQ(mbf.nsetbit(), alpha_occ.size() + beta_occ.size());

    uintv_t chk_inds;

    auto& occs = mbf.m_decoded.m_spin_sym_occs.get();
    auto& vacs = mbf.m_decoded.m_spin_sym_vacs.get();
    /*
     * alpha occupied
     */
    ASSERT_EQ(occs.size({0, 0}), 3);
    chk_inds = occs[{0, 0}];
    ASSERT_EQ(chk_inds, uintv_t({0, 4, 9}));

    ASSERT_EQ(occs.size({0, 1}), 2);
    chk_inds = occs[{0, 1}];
    ASSERT_EQ(chk_inds, uintv_t({1, 2}));

    ASSERT_EQ(occs.size({0, 2}), 0);
    chk_inds = occs[{0, 2}];
    ASSERT_EQ(chk_inds, uintv_t({}));

    ASSERT_EQ(occs.size({0, 3}), 1);
    chk_inds = occs[{0, 3}];
    ASSERT_EQ(chk_inds, uintv_t({7}));

    /*
     * alpha vacant
     */
    ASSERT_EQ(vacs.size({0, 0}), 1);
    chk_inds = vacs[{0, 0}];
    ASSERT_EQ(chk_inds, uintv_t({5}));

    ASSERT_EQ(vacs.size({0, 1}), 2);
    chk_inds = vacs[{0, 1}];
    ASSERT_EQ(chk_inds, uintv_t({6, 8}));

    ASSERT_EQ(vacs.size({0, 2}), 0);
    chk_inds = vacs[{0, 2}];
    ASSERT_EQ(chk_inds, uintv_t({}));

    ASSERT_EQ(vacs.size({0, 3}), 1);
    chk_inds = vacs[{0, 3}];
    ASSERT_EQ(chk_inds, uintv_t({3}));

    /*
     * beta occupied
     * beta_occ{2, 4, 5, 6, 7, 8} -> 12, 14, 15, 16, 17, 18
     */
    ASSERT_EQ(occs.size({1, 0}), 2);
    chk_inds = occs[{1, 0}];
    ASSERT_EQ(chk_inds, uintv_t({14, 15}));

    ASSERT_EQ(occs.size({1, 1}), 3);
    chk_inds = occs[{1, 1}];
    ASSERT_EQ(chk_inds, uintv_t({12, 16, 18}));

    ASSERT_EQ(occs.size({1, 2}), 0);
    chk_inds = occs[{1, 2}];
    ASSERT_EQ(chk_inds, uintv_t({}));

    ASSERT_EQ(occs.size({1, 3}), 1);
    chk_inds = occs[{1, 3}];
    ASSERT_EQ(chk_inds, uintv_t({17}));

    /*
     * beta vacant
     */
    ASSERT_EQ(vacs.size({1, 0}), 2);
    chk_inds = vacs[{1, 0}];
    ASSERT_EQ(chk_inds, uintv_t({10, 19}));

    ASSERT_EQ(vacs.size({1, 1}), 1);
    chk_inds = vacs[{1, 1}];
    ASSERT_EQ(chk_inds, uintv_t({11}));

    ASSERT_EQ(vacs.size({1, 2}), 0);
    chk_inds = vacs[{1, 2}];
    ASSERT_EQ(chk_inds, uintv_t({}));

    ASSERT_EQ(vacs.size({1, 3}), 1);
    chk_inds = vacs[{1, 3}];
    ASSERT_EQ(chk_inds, uintv_t({13}));

    /*
     * all spinorb labels except those correcsponing to irrep "Y" have nonempty pairs
     */
    auto& nonempty_pair_labels = mbf.m_decoded.m_nonempty_pair_labels.get();
    chk_inds = {0, 1, 3, 4, 5, 7};
    ASSERT_EQ(chk_inds, nonempty_pair_labels);
}

TEST(DecodedMbf, Bosons) {
    const uint_t nmode = 20;
    buffered::BosOnv mbf(nmode);
    mbf = {0, 1, 0, 3, 1, 0, 1, 2, 0, 0, 1, 3, 0, 0, 0, 4, 1, 0, 2, 1};
    uintv_t chk_inds;

    chk_inds = {1, 3, 3, 3, 4, 6, 7, 7, 10, 11, 11, 11, 15, 15, 15, 15, 16, 18, 18, 19};
    ASSERT_EQ(mbf.m_decoded.m_expanded.get(), chk_inds);

    chk_inds = {1, 3, 4, 6, 7, 10, 11, 15, 16, 18, 19};
    ASSERT_EQ(mbf.m_decoded.m_occ_modes.get(), chk_inds);
}

TEST(DecodedMbf, Holstein) {
    const uint_t nsite = 15;
    buffered::FrmBosOnv mbf(nsite, nsite);
    mbf.m_frm = {{0, 1, 5, 7, 10, 11, 14}, {0, 3, 6, 7, 11, 12, 14}};
    mbf.m_bos = {0, 1, 2, 3, 1, 0, 1, 2, 0, 1, 3, 0, 0, 1, 2};

    uintv_t chk_inds = {1, 3, 6, 7, 10, 14};
    ASSERT_EQ(mbf.m_decoded.m_occ_sites_nonzero_bosons.get(), chk_inds);
}
