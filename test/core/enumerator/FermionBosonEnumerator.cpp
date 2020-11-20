//
// Created by rja on 05/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/enumerator/FermiBosOnvEnumerator.h"

TEST(FermionBosonEnumerator, SpinNonCon){
    const size_t nsite=4, nelec=4, nmode=4, occ_cutoff=3;
    FermiBosOnvEnumerator enumerator(nsite, nelec, nmode, occ_cutoff);
    size_t i = ~0ul;
    elements::FermiBosOnv config(nsite, nmode);
    while(enumerator.next(config, i)){}
    ASSERT_EQ(i, ci_utils::fermion_dim(nsite, nelec)*ci_utils::boson_dim(nmode, occ_cutoff));
}

TEST(FermionBosonEnumerator, SpinCon){
    const size_t nsite=4, nelec=4, nmode=4, occ_cutoff=3;
    int spin=0;
    FermiBosOnvEnumerator enumerator(nsite, nelec, spin, nmode, occ_cutoff);
    size_t i = ~0ul;
    elements::FermiBosOnv config(nsite, nmode);
    while(enumerator.next(config, i)){}
    ASSERT_EQ(i, ci_utils::fermion_dim(nsite, nelec, spin)*ci_utils::boson_dim(nmode, occ_cutoff));
}