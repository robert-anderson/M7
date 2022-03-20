//
// Created by Robert John Anderson on 2020-01-18.
//

#include <gtest/gtest.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/field/Mbf.h>
#include "M7_lib/caches/DecodedDeterminants.h"

#ifdef ENABLE_COMPLEX
TEST(FermionHamiltonian, DhfEnergy) {
    const auto benchmark = -14.354220448530139;
    FermionHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_nsite);
    onv = {0, 1, ham.m_nsite, ham.m_nsite + 1};
    auto elem = ham.get_element(onv);
    ASSERT_TRUE(consts::floats_equal(consts::real(elem), benchmark));
    ASSERT_TRUE(consts::float_nearly_zero(consts::imag(elem), 1e-14));
    ASSERT_TRUE(consts::floats_equal(ham.get_energy(onv), benchmark));
}

TEST(FermionHamiltonian, DhfBrillouinTheorem) {
    FermionHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv hf_det(ham.m_nsite);
    hf_det = {0, 1, ham.m_nsite, ham.m_nsite + 1};

    OccOrbs occs(hf_det);
    VacOrbs vacs(hf_det);

    buffered::FrmOnv excited(ham.m_nsite);

    conn::FrmOnv conn(hf_det.m_nsite);

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[iocc];
            conn.set(occ, vac);
            ASSERT_TRUE(consts::float_is_zero(ham.get_element(hf_det, conn)));
        }
    }
}
#endif

TEST(FermionHamiltonian, RhfEnergy) {
    const auto benchmark = -108.76171800006861;
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_N2_6o6e/FCIDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    defs::inds chk_orbsyms = {0, 2, 1, 5, 6, 4};
    ASSERT_EQ(ham.m_frm->m_point_group_map.m_site_irreps, chk_orbsyms);
    ASSERT_TRUE(ham.m_frm->m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_bd.m_nsite);
    mbf::set_aufbau_mbf(onv, ham);
    auto elem = ham.get_element(onv);
    ASSERT_FLOAT_EQ(consts::real(elem), benchmark);
    ASSERT_FLOAT_EQ(consts::imag(elem), 1e-14);
    ASSERT_FLOAT_EQ(ham.get_energy(onv), benchmark);
}

TEST(FermionHamiltonian, RhfBrillouinTheorem) {
    fciqmc_config::Document opts;
    opts.m_hamiltonian.m_fermion.m_fcidump.m_path = defs::assets_root + "/RHF_N2_6o6e/FCIDUMP";
    opts.verify();
    Hamiltonian ham(opts.m_hamiltonian);
    ASSERT_TRUE(ham.m_frm->m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_bd.m_nsite);
    mbf::set_aufbau_mbf(onv, ham);

    OccOrbs occs(onv);
    VacOrbs vacs(onv);

    buffered::FrmOnv excited(ham.m_bd.m_nsite);

    conn::FrmOnv conn(onv.m_nsite);

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[iocc];
            conn.clear();
            conn.add(occ, vac);
            ASSERT_FLOAT_EQ(ham.m_frm->get_element_1100(onv, conn), 0.0);
        }
    }
}