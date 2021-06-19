//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/table/BufferedFields.h>
#include "src/core/basis/DecodedDeterminant.h"
#include "gtest/gtest.h"

TEST(DecodedDeterminant, CopyAndMove){
    const size_t n = 6;
    OccupiedOrbitals occorbs(n);
    ASSERT_EQ(occorbs.inds().capacity(), 2*n);
    auto occorbs_cpy = occorbs;
    ASSERT_EQ(occorbs_cpy.inds().capacity(), 2*n);
    auto occorbs_mov = std::move(occorbs);
    ASSERT_EQ(occorbs_mov.inds().capacity(), 2*n);
}

TEST(DecodedDeterminant, Occupation){
    buffered::Onv<0> onv(50);
    defs::inds occ{0, 1, 4, 7, 32, 50, 51, 54, 60, 89, 99};
    onv = occ;
    defs::inds vac;
    auto iter = occ.begin();
    for (size_t i=0ul; i < onv.nbit(); ++i){
        if (iter!=occ.end() && i==*iter) iter++;
        else vac.push_back(i);
    }

    OccupiedOrbitals occorbs(onv);
    ASSERT_TRUE(std::equal(occ.begin(), occ.end(), occorbs.inds().begin()));

    VacantOrbitals vacorbs(onv);
    ASSERT_TRUE(std::equal(vac.begin(), vac.end(), vacorbs.inds().begin()));
}

TEST(DecodedDeterminant, SymmDecoded){

    buffered::Onv<0> onv(10);
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

    SymmOccupiedOrbitals occorbs(grp_map);
    SymmVacantOrbitals vacorbs(grp_map);

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