//
// Created by rja on 05/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/enumerator/FrmBosOnvEnumerator.h"

TEST(FermionBosonEnumerator, SpinNonCon){
    const size_t nsite=4, nelec=4, nmode=4, occ_cutoff=3;
    const BasisDims bd = {nsite, nmode};
    FrmBosOnvEnumerator enumerator(bd, nelec, occ_cutoff);
    size_t i = ~0ul;
    buffered::FrmBosOnv mbf(bd);
    while(enumerator.next(mbf, i)){}
    ASSERT_EQ(i, ci_utils::fermion_dim(nsite, nelec)*ci_utils::boson_dim(nmode, occ_cutoff, false));
}

TEST(FermionBosonEnumerator, SpinCon){
    const size_t nsite=4, nelec=4, nmode=4, occ_cutoff=3;
    const BasisDims bd = {nsite, nmode};
    int spin=0;
    FrmBosOnvEnumerator enumerator(bd, nelec, spin, occ_cutoff);
    size_t i = ~0ul;
    buffered::FrmBosOnv mbf(bd);
    while(enumerator.next(mbf, i)){}
    ASSERT_EQ(i, ci_utils::fermion_dim(nsite, nelec, spin)*ci_utils::boson_dim(nmode, occ_cutoff, false));
}