//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/FciInitializer.h>
#include "test_core/defs.h"

TEST(FciInitializer, N2) {
    //GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"}, true);
    HubbardFrmHam frm_ham(10.0, lattice::make("ortho", {8}, {-1}));
    Hamiltonian ham(&frm_ham);
    FciInitializer init(ham, -25.2245);
    ASSERT_NEARLY_EQ(init.m_eval, -5.834322635772544);
}