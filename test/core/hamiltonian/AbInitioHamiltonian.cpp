//
// Created by Robert John Anderson on 2020-01-18.
//

#include <gtest/gtest.h>
#include <src/core/fermion/DecodedDeterminant.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"

TEST(AbInitioHamiltonian, DhfEnergy){
    AbInitioHamiltonian ham(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
    Determinant hf_det(ham.nsite());
    hf_det.set(defs::inds{0, 1, ham.nsite(), ham.nsite()+1});
    auto elem = ham.get_element_0(hf_det);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), -14.354220448530139));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
}

TEST(AbInitioHamiltonian, DhfBrillouinTheorem) {
    AbInitioHamiltonian ham(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
    Determinant hf_det(ham.nsite());
    hf_det.set(defs::inds{0, 1, ham.nsite(), ham.nsite()+1});
    size_t removed, inserted;

    OccupiedOrbitals occs(hf_det);
    VacantOrbitals vacs(hf_det);

    Determinant excited(ham.nsite());

    for (size_t iocc=0ul; iocc<occs.m_nind; ++iocc){
        const auto &occ = occs.m_inds[iocc];
        for (size_t ivac=0ul; ivac<vacs.m_nind; ++ivac) {
            const auto &vac = vacs.m_inds[iocc];
            ASSERT_TRUE(consts::float_is_zero(ham.get_element_1(hf_det, occ, vac)));
            excited = hf_det;
            excited.excite(vac, occ);
            AntisymConnection connection(hf_det, excited);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element(hf_det, connection)));
            connection.update(excited, hf_det);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element(excited, connection)));
        }
    }
}