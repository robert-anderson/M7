//
// Created by anderson on 10/06/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"

TEST(GeneralFrmBosHam, Element1101hc) {
    GeneralLadderHam frmbos_ham({PROJECT_ROOT"/assets/Hubbard_U4_3site/EBDUMP_GENERAL"}, true, 1);
    const auto& basis = frmbos_ham.m_basis.m_frm;
    // 0.29600971859745534    3    1    2
    const auto coeff = 0.29600971859745534;
    ASSERT_EQ(frmbos_ham.get_coeff_1110(2, basis.ispinorb(0, 0), basis.ispinorb(0, 1)), coeff);
    // check other spin channel:
    ASSERT_EQ(frmbos_ham.get_coeff_1110(2, basis.ispinorb(1, 0), basis.ispinorb(1, 1)), coeff);
    // and conjugate:
    ASSERT_EQ(frmbos_ham.get_coeff_1101(2, basis.ispinorb(0, 1), basis.ispinorb(0, 0)), coeff);
    ASSERT_EQ(frmbos_ham.get_coeff_1101(2, basis.ispinorb(1, 1), basis.ispinorb(1, 0)), coeff);
}