//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/table/BufferedFields.h>
#include "src/core/basis/DecodedDeterminant.h"
#include "gtest/gtest.h"

TEST(DecodedDeterminant, Occupation){
    buffered::FermionOnv det(50);
    defs::inds occ{0, 1, 4, 7, 32, 50, 51, 54, 60, 89, 99};
    det = occ;
    defs::inds vac;
    auto iter = occ.begin();
    for (size_t i=0ul; i<det.m_nbit; ++i){
        if (iter!=occ.end() && i==*iter) iter++;
        else vac.push_back(i);
    }

    OccupiedOrbitals occorbs(det);
    ASSERT_TRUE(std::equal(occ.begin(), occ.end(), occorbs.inds().begin()));

    VacantOrbitals vacorbs(det);
    ASSERT_TRUE(std::equal(vac.begin(), vac.end(), vacorbs.inds().begin()));
}