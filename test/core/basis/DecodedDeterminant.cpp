//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/foreach/ForeachVirtual.h>
#include "src/core/caches/CachedOrbs.h"
#include "gtest/gtest.h"

TEST(DecodedDeterminant, SimpleOccAndVac){
    using namespace decoded_mbf::spinorbs;
    buffered::FrmOnv mbf({50, 0});
    defs::inds occ{0, 1, 4, 7, 32, 50, 51, 54, 60, 89, 99};
    mbf = occ;
    defs::inds vac;
    auto iter = occ.begin();
    for (size_t i=0ul; i < mbf.nbit(); ++i){
        if (iter!=occ.end() && i==*iter) iter++;
        else vac.push_back(i);
    }

    SimpleOccs occorbs;
    occorbs.update(mbf);
    ASSERT_TRUE(std::equal(occ.begin(), occ.end(), occorbs.inds().begin()));

    SimpleVacs vacorbs;
    vacorbs.update(mbf);
    ASSERT_TRUE(std::equal(vac.begin(), vac.end(), vacorbs.inds().begin()));

    /*
     * for a small number of occupied orbs, run through all possible arrangements
     */
    const size_t noccorb = 3;
    auto occ_fn = [&mbf, &occorbs](size_t iiter, const defs::inds &value) {
        mbf.zero();
        mbf = value;
        occorbs.update(mbf);
        ASSERT_EQ(occorbs.inds(), value);
    };
    foreach_virtual::rtnd::lambda::Ordered<> occ_foreach(occ_fn, mbf.m_nspinorb, noccorb);
    occ_foreach.loop();

    /*
     * for a small number of vacant orbs, run through all possible arrangements
     */
    const size_t nvacorb = 3;
    auto vac_fn = [&mbf, &vacorbs](size_t iiter, const defs::inds &value) {
        mbf.set();
        for (auto i: value) mbf.clr(i);
        vacorbs.update(mbf);
        ASSERT_EQ(vacorbs.inds(), value);
    };
    foreach_virtual::rtnd::lambda::Ordered<> vac_foreach(vac_fn, mbf.m_nspinorb, nvacorb);
    vac_foreach.loop();
}

TEST(DecodedDeterminant, SymmDecoded){
    using namespace decoded_mbf::spinorbs;

    buffered::FrmOnv onv(10);
    defs::inds alpha_occ{0, 1, 2, 4, 7, 9};
    defs::inds beta_occ{2, 4, 5, 6, 7, 8};
    onv = {alpha_occ, beta_occ};
    ASSERT_EQ(onv.nsetbit(), alpha_occ.size()+beta_occ.size());

    // arbitrary, fictitious group
    AbelianGroup grp({"X", "Y", "Z"}, [](const size_t& iirrep, const size_t& jirrep){
        return (iirrep+jirrep)%3;
    });

    /*                                         0  1  2  3  4  5  6  7  8  9
     * occupied orbitals (alpha):              o  o  o  /  o  /  /  o  /  o
     * occupied orbitals (beta ):              /  /  o  /  o  o  o  o  o  /
     */
    AbelianGroupMap grp_map(grp, {0, 1, 1, 2, 0, 0, 1, 2, 1, 0});
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

    ASSERT_EQ(grp_map.m_nsite, 10);

    SpinSymOccs occorbs(grp_map);
    SpinSymVacs vacorbs(grp_map);

    occorbs.update(onv);
    vacorbs.update(onv);

    defs::inds chk_inds;

    /*
     * alpha occupied
     */
    ASSERT_EQ(occorbs.size(0), 3);
    chk_inds = occorbs[{0, 0}];
    ASSERT_EQ(chk_inds, defs::inds({0, 4, 9}));
    ASSERT_EQ(occorbs.size(1), 2);
    chk_inds = occorbs[{0, 1}];
    ASSERT_EQ(chk_inds, defs::inds({1, 2}));
    ASSERT_EQ(occorbs.size(2), 1);
    chk_inds = occorbs[{0, 2}];
    ASSERT_EQ(chk_inds, defs::inds({7}));

    /*
     * alpha vacant
     */
    ASSERT_EQ(vacorbs.size(0), 1);
    chk_inds = vacorbs[{0, 0}];
    ASSERT_EQ(chk_inds, defs::inds({5}));
    ASSERT_EQ(vacorbs.size(1), 2);
    chk_inds = vacorbs[{0, 1}];
    ASSERT_EQ(chk_inds, defs::inds({6, 8}));
    ASSERT_EQ(vacorbs.size(2), 1);
    chk_inds = vacorbs[{0, 2}];
    ASSERT_EQ(chk_inds, defs::inds({3}));


    /*
     * beta occupied
     * beta_occ{2, 4, 5, 6, 7, 8} -> 12, 14, 15, 16, 17, 18
     */
    ASSERT_EQ(occorbs.size(3), 2);
    chk_inds = occorbs[{1, 0}];
    ASSERT_EQ(chk_inds, defs::inds({14, 15}));
    ASSERT_EQ(occorbs.size(4), 3);
    chk_inds = occorbs[{1, 1}];
    ASSERT_EQ(chk_inds, defs::inds({12, 16, 18}));
    ASSERT_EQ(occorbs.size(5), 1);
    chk_inds = occorbs[{1, 2}];
    ASSERT_EQ(chk_inds, defs::inds({17}));

    /*
     * beta vacant
     */
    ASSERT_EQ(vacorbs.size(3), 2);
    chk_inds = vacorbs[{1, 0}];
    ASSERT_EQ(chk_inds, defs::inds({10, 19}));
    ASSERT_EQ(vacorbs.size(4), 1);
    chk_inds = vacorbs[{1, 1}];
    ASSERT_EQ(chk_inds, defs::inds({11}));
    ASSERT_EQ(vacorbs.size(5), 1);
    chk_inds = vacorbs[{1, 2}];
    ASSERT_EQ(chk_inds, defs::inds({13}));
}