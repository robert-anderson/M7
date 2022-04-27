//
// Created by Robert J. Anderson on 08/03/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"
#include "M7_lib/hamiltonian/Hamiltonian.h"

TEST(InteractingBosGasBosHam, DiagonalMatrixElements) {
    const size_t nwave = 3;
    conf::Hamiltonian opts(nullptr);
    opts.m_boson.m_interacting_bose_gas.m_ndim = 1;
    opts.m_boson.m_interacting_bose_gas.m_nwave = nwave;
    opts.m_boson.m_interacting_bose_gas.m_ek_scale = 1.0;
    opts.m_boson.m_nboson = 3;
    Hamiltonian ham_src(opts);
    buffered::BosOnv mbf(ham_src.m_hs);
    mbf = {0, 0, 0, 0, 0, 0, 3};
    defs::ham_t helem;
    helem = ham_src.get_element(mbf);
    ASSERT_TRUE(consts::nearly_equal(helem, 5.07));
}