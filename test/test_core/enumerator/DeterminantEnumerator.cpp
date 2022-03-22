//
// Created by jhalson on 21/10/2020.
//

#include <gtest/gtest.h>
#include "M7_lib/enumerator/FermionOnvEnumerator.h"
#include "M7_lib/table/BufferedFields.h"


TEST(SpinNonConDetEnumerator, SimpleConstruction){
    size_t nsite = 4, nelec = 4, idet=~0ul;
    FermionOnvEnumerator enumerator(nsite, nelec);
    buffered::FrmOnv det(nsite);
    while(enumerator.next(det, idet)){}
    ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec));
}

TEST(SpinConDetEnumerator, EnumerateSpinZero){
    size_t nsite = 4, nelec = 4, idet=~0ul;
    int ms2 = 0;
    FermionOnvEnumerator enumerator(nsite, nelec, ms2);
    buffered::FrmOnv det(nsite);
    while(enumerator.next(det, idet)){
        ASSERT_EQ(det.ms2(), ms2);
    }
    ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec, ms2));
}


TEST(SpinConDetEnumerator, EnumerateSpinOdd){
    size_t nsite = 5, nelec = 5;
    for(auto ms2 : {-5, -3, -1, 1, 3, 5}){
        size_t idet=~0ul;
        FermionOnvEnumerator scde(nsite, nelec, ms2);
        buffered::FrmOnv det(nsite);
        while(scde.next(det, idet)){
            ASSERT_EQ(det.ms2(), ms2);
        }
        ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec, ms2));
    }
}

TEST(SpinConDetEnumerator, EnumerateSpinEven){
    size_t nsite = 6, nelec = 6;
    for(auto ms2 : {-6, -4, -2, 0, 2, 4, 6}){
        size_t idet=~0ul;
        FermionOnvEnumerator scde(nsite, nelec, ms2);
        buffered::FrmOnv det(nsite);
        while(scde.next(det, idet)){
            ASSERT_EQ(det.ms2(), ms2);
        }
        ASSERT_EQ(idet, ci_utils::fermion_dim(nsite, nelec, ms2));
    }
}