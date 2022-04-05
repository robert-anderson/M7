//
// Created by jhalson on 30/09/2020.
//

#include <M7_lib/table/BufferedFields.h>
#include "gtest/gtest.h"
#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"
#include "M7_lib/hamiltonian/bos/BosHam.h"

#if 0
TEST(BosonCouplings, Element_b0) {
    size_t nboson_cutoff = 4, nsite = 4;
    defs::ham_t v = 0.5, omega = 0.025;
    conn::FrmBosOnv conn(nsite);

    BosonCouplings bc(nsite, nboson_cutoff, v);
    BosonHamiltonian bos_ham(nsite, nboson_cutoff, omega);

    buffered::FrmBosOnv src(nsite);
    ASSERT_TRUE(src.m_row);
    ASSERT_TRUE(src.m_frm.belongs_to_row());
    ASSERT_TRUE(src.m_frm.m_row);
    ASSERT_TRUE(src.m_bos.belongs_to_row());
    ASSERT_TRUE(src.m_bos.m_row);
    ASSERT_EQ(src.m_row, src.m_frm.m_row);
    ASSERT_EQ(src.m_row, src.m_bos.m_row);
    buffered::FrmBosOnv dst(nsite);

    src = {{1, 2, 3, 4},
           {0, 0, 0, 0}};
    dst = {{1, 2, 3, 4},
           {0, 0, 0, 0}};
    conn.connect(src, dst);
    defs::ham_t helement;

    helement = bc.get_element(src, conn);
    ASSERT_EQ(helement, 0.0);
    helement = bos_ham.get_element(src.m_bos, conn.m_bos);
    ASSERT_EQ(helement, 0.0);

    src = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    dst = {{1, 2, 3, 4},
           {2, 4, 0, 1}};

    conn.connect(src, dst);

    helement = bc.get_element(src, conn);
    // no change in boson occupation, so no coupling
    ASSERT_EQ(helement, 0.0);
    helement = bos_ham.get_element(src.m_bos, conn.m_bos);
    ASSERT_EQ(helement, 7.0 * omega);
}

TEST(BosonCouplings, Element_f0_b1){
    size_t nboson_cutoff = 4, nsite = 4;
    defs::ham_t v = 0.5, omega = 0.025;
    conn::FrmBosOnv conn(nsite);

    BosonCouplings bc(nsite, nboson_cutoff, v);
    BosonHamiltonian bos_ham(nsite, nboson_cutoff, omega);

    buffered::FrmBosOnv src(nsite);
    buffered::FrmBosOnv dst(nsite);

    src = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    dst = {{1, 2, 3, 4},
           {2, 4, 0, 2}};
    conn.connect(src, dst);

    ASSERT_EQ(conn.m_bos.size(), 1ul);
    ASSERT_EQ(conn.m_bos.m_cre[0].m_imode, 3ul);
    ASSERT_EQ(conn.m_bos.m_cre[0].m_nop, 1ul);
    ASSERT_EQ(conn.m_bos.m_ann.size(), 0ul);

    defs::ham_t helement;

    helement = bc.get_element(src, conn);
    ASSERT_EQ(helement, v*std::sqrt(2.0));
    helement = bos_ham.get_element(src.m_bos, conn.m_bos);
    ASSERT_EQ(helement, 0.0);

    src = {{1, 2, 3, 4},
           {2, 6, 0, 1}};
    dst = {{1, 2, 3, 4},
           {2, 5, 0, 1}};
    conn.connect(src, dst);

    ASSERT_EQ(conn.m_bos.size(), 1ul);
    ASSERT_EQ(conn.m_bos.m_ann[0].m_imode, 1ul);
    ASSERT_EQ(conn.m_bos.m_ann[0].m_nop, 1ul);
    ASSERT_EQ(conn.m_bos.m_cre.size(), 0ul);

    helement = bc.get_element(src, conn);
    ASSERT_EQ(helement, v*std::sqrt(6.0));
    helement = bos_ham.get_element(src.m_bos, conn.m_bos);
    ASSERT_EQ(helement, 0.0);


    src = {{1, 2, 3, 4},
           {2, 6, 0, 1}};
    dst = {{1, 2, 3, 4},
           {2, 7, 0, 1}};
    conn.connect(src, dst);

    ASSERT_EQ(conn.m_bos.size(), 1ul);
    ASSERT_EQ(conn.m_bos.m_cre[0].m_imode, 1ul);
    ASSERT_EQ(conn.m_bos.m_cre[0].m_nop, 1ul);
    ASSERT_EQ(conn.m_bos.m_ann.size(), 0ul);

    helement = bc.get_element(src, conn);
    ASSERT_EQ(helement, v*std::sqrt(7.0));
    helement = bos_ham.get_element(src.m_bos, conn.m_bos);
    ASSERT_EQ(helement, 0.0);
}

TEST(BosonCouplings, Element_f1_b1){
    size_t nboson_cutoff = 4, nsite = 4;
    defs::ham_t v = 0.5, omega = 0.025;
    conn::FrmBosOnv conn(nsite);

    BosonCouplings bc(nsite, nboson_cutoff, v);
    BosonHamiltonian bos_ham(nsite, nboson_cutoff, omega);

    buffered::FrmBosOnv src(nsite);
    buffered::FrmBosOnv dst(nsite);

    src = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    dst = {{1, 2, 3, 5},
           {2, 4, 0, 2}};

    conn.connect(src, dst);
    defs::ham_t helement;
    helement = bc.get_element(src, conn);
    ASSERT_EQ(helement, 0.0);
}
#endif