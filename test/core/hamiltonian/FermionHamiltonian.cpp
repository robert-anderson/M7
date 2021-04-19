//
// Created by Robert John Anderson on 2020-01-18.
//

#include <gtest/gtest.h>
#include "src/core/basis/DecodedDeterminant.h"
#include "src/core/hamiltonian/FermionHamiltonian.h"

#ifdef ENABLE_COMPLEX
TEST(FermionHamiltonian, DhfEnergy) {
    const auto benchmark = -14.354220448530139;
    FermionHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.spin_conserving());
    buffered::FermionOnv fonv(ham.nsite());
    fonv = {0, 1, ham.nsite(), ham.nsite() + 1};
    auto elem = ham.get_element_0(fonv);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), benchmark));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
    ASSERT_TRUE(consts::floats_equal(ham.get_energy(fonv), benchmark));
}

TEST(FermionHamiltonian, DhfBrillouinTheorem) {
    FermionHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.spin_conserving());
    buffered::FermionOnv hf_det(ham.nsite());
    hf_det = {0, 1, ham.nsite(), ham.nsite() + 1};

    OccupiedOrbitals occs(hf_det);
    VacantOrbitals vacs(hf_det);

    buffered::FermionOnv excited(ham.nsite());

    AntisymFermionOnvConnection connection(hf_det);

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[iocc];
            connection.zero();
            connection.add(occ, vac);
            connection.apply(hf_det);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element_1(connection)));
            /*
            excited = hf_det;
            excited.excite(vac, occ);
            AntisymFermionOnvConnection connection(hf_det, excited);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element(hf_det, connection)));
            connection.connect(excited, hf_det);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element(excited, connection)));
             */
        }
    }
}
#endif

TEST(FermionHamiltonian, RhfEnergy) {
    const auto benchmark = -108.65146156994338;
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv fonv(ham.nsite());
    for (size_t i=0ul; i<ham.nelec()/2; ++i){fonv.set({0, i}); fonv.set({1, i});}

    auto elem = ham.get_element_0(fonv);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), benchmark));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
    ASSERT_TRUE(consts::floats_equal(ham.get_energy(fonv), benchmark));
}

TEST(FermionHamiltonian, RhfBrillouinTheorem) {
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    buffered::FermionOnv fonv(ham.nsite());
    fonv = {0, 1, 2,  6, 7, 8};

    OccupiedOrbitals occs(fonv);
    VacantOrbitals vacs(fonv);

    buffered::FermionOnv excited(ham.nsite());

    AntisymFermionOnvConnection connection(fonv);

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[iocc];
            connection.zero();
            connection.add(occ, vac);
            connection.apply(fonv);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element_1(connection)));
        }
    }
}