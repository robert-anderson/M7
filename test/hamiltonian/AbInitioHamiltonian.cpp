//
// Created by Robert John Anderson on 2020-01-18.
//


#include <gtest/gtest.h>
#include "src/enumerators/BitfieldEnumerator.h"
#include "src/hamiltonian/AbInitioHamiltonian.h"

TEST(AbInitioHamiltonian, DhfEnergy){
    AbInitioHamiltonian ham(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
    Determinant hf_det(ham.nsite());
    hf_det.m_bitfields[0].set(defs::inds{0,1});
    hf_det.m_bitfields[1].set(defs::inds{0,1});
    auto elem = ham.get_element_0(hf_det);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), -14.354220448530139));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
}

TEST(AbInitioHamiltonian, DhfBrillouinTheorem) {
    AbInitioHamiltonian ham(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
    Determinant hf_det(ham.nsite());
    hf_det.m_bitfields[0].set(defs::inds{0,1});
    hf_det.m_bitfields[1].set(defs::inds{0,1});
    size_t removed, inserted;
    DeterminantSetEnumerator occupied(hf_det);
    while (occupied.next(removed)){
        {
            DeterminantClrEnumerator unoccupied(hf_det);
            while (unoccupied.next(inserted)){
                ASSERT_TRUE(consts::float_is_zero(ham.get_element_1(hf_det, removed, inserted)));
            }
        }
    }
}
