//
// Created by rja on 07/11/22.
//

#include "test_core/defs.h"
#include "M7_lib/bilinear/FockRdm4.h"

TEST(FockMatrix, LoadingFromMolcasHdf5) {
    FockMatrix fock(5, PROJECT_ROOT"/assets/HF_RDMs/fock.h5");
    ASSERT_NEAR_EQ(fock(0, 0), -6.575519e-01);
    ASSERT_NEAR_EQ(fock(1, 1), -6.575519e-01);
    ASSERT_NEAR_EQ(fock(2, 2), -5.884330e-01);
    ASSERT_NEAR_EQ(fock(2, 3), -6.368367e-03);
    ASSERT_NEAR_EQ(fock(3, 2), -6.368367e-03);
    ASSERT_NEAR_EQ(fock(4, 4),  6.051609e-01);
    ASSERT_NEAR_EQ(fock(3, 4), -6.940660e-03);
    ASSERT_NEAR_EQ(fock(4, 3), -6.940660e-03);
}