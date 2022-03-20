//
// Created by rja on 21/11/2021.
//

#include <src/core/hamiltonian/GeneralBosHam.h>
#include "src/core/config/FciqmcConfig.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "gtest/gtest.h"

TEST(BosonHamiltonian, Coefficients) {
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/LandauLevels_5_5_15/BOSDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    ASSERT_EQ(ham.m_bd.m_nmode, 5ul);
    ASSERT_EQ(ham.nboson(), 5ul);
    //0.2209708691 2 4 5 3
    auto& bos_ham = dynamic_cast<const GeneralBosHam&>(*ham.m_bos);
    ASSERT_FLOAT_EQ(bos_ham.m_coeffs_2.get(1, 3, 4, 2), 0.2209708691);
    ASSERT_FLOAT_EQ(bos_ham.m_coeffs_2.phys_element(1, 4, 3, 2), 0.2209708691);
    //0.1530931089 5 3 1 3
    ASSERT_FLOAT_EQ(bos_ham.m_coeffs_2.get(4, 2, 0, 2), 0.1530931089);
    ASSERT_FLOAT_EQ(bos_ham.m_coeffs_2.phys_element(4, 0, 2, 2), 0.1530931089);

    ASSERT_FLOAT_EQ(bos_ham.m_coeffs_2.get(0, 3, 0, 3), 0.0);
    ASSERT_FLOAT_EQ(bos_ham.m_coeffs_2.phys_element(0, 0, 3, 3), 0.0);
}

TEST(BosonHamiltonian, DiagonalMatrixElements) {
    /*
     *  ONV (L=15)       diag. matrix element
     *  [0. 0. 0. 5. 0.] 3.125
     *  [0. 0. 1. 3. 1.] 4.921875
     *  [0. 0. 2. 1. 2.] 4.8671875
     *  [0. 1. 0. 2. 2.] 4.3984375
     *  [0. 1. 1. 0. 3.] 3.9140625
     *  [1. 0. 0. 1. 3.] 3.0859375
     */
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/LandauLevels_5_5_15/BOSDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    buffered::BosOnv onv(ham.m_bd);
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

TEST(BosonHamiltonian, OffDiagonalMatrixElements) {
    /*
     *  ONVs
     *  [0, 0, 0, 5, 0]
     *  [0, 0, 1, 3, 1]
     *  [0, 0, 2, 1, 2]
     *  [0, 1, 0, 2, 2]
     *  [0, 1, 1, 0, 3]
     *  [1, 0, 0, 1, 3]
     *
     * [[3.125      1.2103073  0.         0.         0.         0.        ]
     *  [1.2103073  4.921875   1.32582521 1.08253175 0.         0.        ]
     *  [0.         1.32582521 4.8671875  0.61237244 1.08253175 0.375     ]
     *  [0.         1.08253175 0.61237244 4.3984375  0.66291261 0.61237244]
     *  [0.         0.         1.08253175 0.66291261 3.9140625  0.4330127 ]
     *  [0.         0.         0.375      0.61237244 0.4330127  3.0859375 ]]
     */
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_boson.m_bosdump.m_path = defs::assets_root + "/LandauLevels_5_5_15/BOSDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    buffered::BosOnv src(ham.m_bd);
    buffered::BosOnv dst(ham.m_bd);
    conn::BosOnv conn(src);

    std::vector<defs::inds> basis =
            {{0, 0, 0, 5, 0},
             {0, 0, 1, 3, 1},
             {0, 0, 2, 1, 2},
             {0, 1, 0, 2, 2},
             {0, 1, 1, 0, 3},
             {1, 0, 0, 1, 3}};

    std::vector<defs::ham_comp_t> h_upper_triangle =
            {1.2103072956898178, 0.0, 0.0, 0.0, 0.0, 1.3258252147247769,
            1.0825317547305484, 0.0, 0.0, 0.6123724356957947,
            1.0825317547305484, 0.37500000000000006, 0.6629126073623882,
            0.6123724356957946, 0.43301270189221935};
    size_t n = 0ul;
    for (size_t i=0ul; i<basis.size(); ++i){
        src = basis[i];
        for (size_t j=i+1; j<basis.size(); ++j) {
            dst = basis[j];
            conn.connect(src, dst);
            ASSERT_FLOAT_EQ(ham.get_element(src, conn), h_upper_triangle[n]);
            // hamiltonian is hermitian
            conn.connect(dst, src);
            ASSERT_FLOAT_EQ(ham.get_element(dst, conn), h_upper_triangle[n]);
            ++n;
        }
    }

    src = basis[0];
    dst = {2, 0, 0, 3, 0};
    // should not be connected due to angular momentum conservation
    conn.connect(src, dst);
    ASSERT_FLOAT_EQ(ham.get_element(src, conn), 0.0);
    conn.connect(dst, src);
    ASSERT_FLOAT_EQ(ham.get_element(dst, conn), 0.0);
}