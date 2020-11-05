//
// Created by jhalson on 21/10/2020.
//

#include <gtest/gtest.h>
#include "src/core/enumerator/DeterminantEnumerator.h"
#include "src/core/field/Elements.h"


TEST(SpinNonConDetEnumerator, SimpleConstruction){
    size_t nsite = 4, nelec = 4, idet=~0ul;
    DeterminantEnumerator enumerator(nsite, nelec);
    elements::Determinant det(nsite);
    while(enumerator.next(det, idet)){}
    ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec));
}

TEST(SpinConDetEnumerator, EnumerateSpinZero){
    size_t nsite = 4, nelec = 4, idet=~0ul;
    int spin = 0;
    DeterminantEnumerator enumerator(nsite, nelec, spin);
    elements::Determinant det(nsite);
    while(enumerator.next(det, idet)){
        ASSERT_EQ(det.spin(), spin);
    }
    ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec, spin));
}


TEST(SpinConDetEnumerator, EnumerateSpinOdd){
    size_t nsite = 5, nelec = 5;
    for(auto spin : {-5, -3, -1, 1, 3, 5}){
        size_t idet=~0ul;
        DeterminantEnumerator scde(nsite, nelec, spin);
        elements::Determinant det(nsite);
        while(scde.next(det, idet)){
            ASSERT_EQ(det.spin(), spin);
        }
        ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec, spin));
    }
}

TEST(SpinConDetEnumerator, EnumerateSpinEven){
    size_t nsite = 6, nelec = 6;
    for(auto spin : {-6, -4,-2, 0, 2, 4, 6}){
        size_t idet=~0ul;
        DeterminantEnumerator scde(nsite, nelec, spin);
        elements::Determinant det(nsite);
        while(scde.next(det, idet)){
            ASSERT_EQ(det.spin(), spin);
        }
        ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec, spin));
    }
}