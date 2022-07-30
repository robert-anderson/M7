//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/wavefunction/FciInitializer.h>
#include <M7_lib/linalg/DenseHamiltonian.h>
#include <M7_lib/hamiltonian/frm/J1J2FrmHam.h>
#include "test_core/defs.h"

TEST(FciInitializer, N2) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP", true});
    Hamiltonian ham(&frm_ham);
    FciInitializer init(ham);
    DenseHamiltonian hmat(ham);
    v_t<ham_t> evals;
    dense::diag(hmat, evals);
    ASSERT_NEARLY_EQ(init.m_eval, evals[0]);
}

TEST(FciInitializer, J1J2) {
    J1J2FrmHam frm_ham(0.25, lattice::make("ortho", {16}, {1}));
    Hamiltonian ham(&frm_ham);
    FciInitializer init(ham, 0.0);
    ASSERT_NEARLY_EQ(init.m_eval, -6.44708);
}