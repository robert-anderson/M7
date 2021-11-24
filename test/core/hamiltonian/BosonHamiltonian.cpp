//
// Created by rja on 21/11/2021.
//

#include "gtest/gtest.h"
#include "src/core/hamiltonian/BosonHamiltonian.h"

TEST(BosonHamiltonian, Coefficients){
    BosonHamiltonian ham(defs::assets_root + "/LandauLevels/BOSDUMP", 12);
    ASSERT_EQ(ham.m_nmode, 5ul);
    ASSERT_EQ(ham.m_nboson, 5ul);
    //0.2209708691 2 4 5 3
    ASSERT_FLOAT_EQ(ham.m_coeffs_2.get(1, 3, 4, 2), 0.2209708691);
    //0.1530931089 5 3 1 3
    ASSERT_FLOAT_EQ(ham.m_coeffs_2.get(4, 2, 0, 2), 0.1530931089);
}

TEST(BosonHamiltonian, DiagonalMatrixElements){
    /*
     *  ONV (L=15)       diag. matrix element
     *  [0. 0. 0. 5. 0.] 3.125
     *  [0. 0. 1. 3. 1.] 4.921875
     *  [0. 0. 2. 1. 2.] 4.8671875
     *  [0. 1. 0. 2. 2.] 4.3984375
     *  [0. 1. 1. 0. 3.] 3.9140625
     *  [1. 0. 0. 1. 3.] 3.0859375
     */
    BosonHamiltonian ham(defs::assets_root + "/LandauLevels/BOSDUMP", 12);
    buffered::BosOnv onv(ham.m_nmode);
    onv = {0, 0, 0, 5, 0};
    ASSERT_FLOAT_EQ(ham.get_element(onv), 3.125);
    onv = {0, 0, 1, 3, 1};
    ASSERT_FLOAT_EQ(ham.get_element(onv), 4.921875);
    onv = {0, 0, 2, 1, 2};
    ASSERT_FLOAT_EQ(ham.get_element(onv), 4.8671875);
    onv = {0, 1, 0, 2, 2};
    ASSERT_FLOAT_EQ(ham.get_element(onv), 4.3984375);
    onv = {0, 1, 1, 0, 3};
    ASSERT_FLOAT_EQ(ham.get_element(onv), 3.9140625);
    onv = {1, 0, 0, 1, 3};
    ASSERT_FLOAT_EQ(ham.get_element(onv), 3.0859375);
}

TEST(BosonHamiltonian, OffDiagonalMatrixElements){
    /*
     *  ONVs
     *  [0. 0. 0. 5. 0.]
     *  [0. 0. 1. 3. 1.]
     *  [0. 0. 2. 1. 2.]
     *  [0. 1. 0. 2. 2.]
     *  [0. 1. 1. 0. 3.]
     *  [1. 0. 0. 1. 3.]
     *
     * [[3.125      1.2103073  0.         0.         0.         0.        ]
     *  [1.2103073  4.921875   1.32582521 1.08253175 0.         0.        ]
     *  [0.         1.32582521 4.8671875  0.61237244 1.08253175 0.375     ]
     *  [0.         1.08253175 0.61237244 4.3984375  0.66291261 0.61237244]
     *  [0.         0.         1.08253175 0.66291261 3.9140625  0.4330127 ]
     *  [0.         0.         0.375      0.61237244 0.4330127  3.0859375 ]]
     */
    BosonHamiltonian ham(defs::assets_root + "/LandauLevels/BOSDUMP", 12);
    buffered::BosOnv src(ham.m_nmode);
    buffered::BosOnv dst(ham.m_nmode);
    src = {0, 0, 0, 5, 0};
    dst = {0, 0, 1, 3, 1};
    conn::BosOnv conn(src);
    conn.connect(src, dst);
    ASSERT_FLOAT_EQ(ham.get_element(src, conn), 1.2103073);
}