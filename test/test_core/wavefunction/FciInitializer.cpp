//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/FciInitializer.h>
#include "test_core/defs.h"

TEST(FciInitializer, N2) {
    //GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"}, true);
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {4, 3}, {1, 1}))                                                                      ;
    Hamiltonian ham(&frm_ham);
    FciInitializer init(ham);
    //ASSERT_NEARLY_EQ(init.m_eval, -108.916561245585);
}