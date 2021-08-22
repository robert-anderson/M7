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
    buffered::FrmOnv onv(ham.m_nsite);
    onv = {0, 1, ham.m_nsite, ham.m_nsite + 1};
    auto elem = ham.get_element(onv);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), benchmark));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
    ASSERT_TRUE(consts::floats_equal(ham.get_energy(onv), benchmark));
}

TEST(FermionHamiltonian, DhfBrillouinTheorem) {
    FermionHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.spin_conserving());
    buffered::FrmOnv hf_det(ham.m_nsite);
    hf_det = {0, 1, ham.m_nsite, ham.m_nsite + 1};

    OccupiedOrbitals occs(hf_det);
    VacantOrbitals vacs(hf_det);

    buffered::FrmOnv excited(ham.m_nsite);

    conn::FrmOnv conn(hf_det.m_nsite);

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[iocc];
            conn.set(occ, vac);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element(hf_det, conn)));
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
    const auto benchmark = -108.76171800006861;
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv fonv(ham.m_nsite);
    for (size_t i=0ul; i<ham.m_nelec/2; ++i){fonv.set({0, i}); fonv.set({1, i});}

    auto elem = ham.get_element(fonv);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), benchmark));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
    ASSERT_TRUE(consts::floats_equal(ham.get_energy(fonv), benchmark));
}

TEST(FermionHamiltonian, RhfBrillouinTheorem) {
    FermionHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv fonv(ham.m_nsite);
    fonv = {0, 1, 2,  6, 7, 8};

    OccupiedOrbitals occs(fonv);
    VacantOrbitals vacs(fonv);

    buffered::FrmOnv excited(ham.m_nsite);

    conn::FrmOnv conn(fonv.m_nsite);

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[iocc];
            conn.clear();
            conn.add(occ, vac);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element_1100(fonv, conn)));
        }
    }
}