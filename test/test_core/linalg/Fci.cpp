//
// Created by rja on 14/07/22.
//

#include "test_core/defs.h"
#include "M7_lib/linalg/Fci.h"

TEST(Fci, Foreach) {
    HubbardFrmHam frm_ham(4.0, lattice::make("ortho", {2, 3}, {0, 0}));
    Hamiltonian ham(&frm_ham);
    auto iters = fci::BasisIters::make(ham);
    buffered::FrmOnv mbf(ham.m_basis);
    auto fn = [&mbf]() {
        std::cout << mbf << std::endl;
    };
    iters.m_single->loop(mbf, fn);
}