//
// Created by Robert John Anderson on 2020-01-18.
//

#include <test_core/defs.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/field/Mbf.h>
#include "M7_lib/caches/DecodedDeterminants.h"

#ifdef ENABLE_COMPLEX_HAM
TEST(FermionHamiltonian, DhfEnergy) {
    const auto benchmark = -14.354220448530139;
    GeneralFrmHam ham({PROJECT_ROOT"/assets/DHF_Be_STO-3G/FCIDUMP"});
    ASSERT_FALSE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_basis);
    onv = {{0, 1}, {0, 1}};
    auto elem = ham.get_element(onv);
    ASSERT_NEAR_EQ(arith::real(elem), benchmark);
    ASSERT_NEAR_ZERO(arith::imag(elem));
    ASSERT_NEAR_EQ(ham.get_energy(onv), benchmark);
}

TEST(FermionHamiltonian, DhfBrillouinTheorem) {
    GeneralFrmHam ham({PROJECT_ROOT"/assets/DHF_Be_STO-3G/FCIDUMP"});
    ASSERT_FALSE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv hf_det(ham.m_basis);
    hf_det = {{0, 1}, {0, 1}};

    buffered::FrmOnv excited(ham.m_basis);
    conn::FrmOnv conn(hf_det.m_basis);

    const auto& occs = hf_det.m_decoded.m_simple_occs.get();
    const auto& vacs = hf_det.m_decoded.m_simple_vacs.get();
    for (uint_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto &occ = occs[iocc];
        for (uint_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            const auto &vac = vacs[ivac];
            conn.m_ann.set(occ);
            conn.m_cre.set(vac);
            ASSERT_NEAR_ZERO(ham.get_element(hf_det, conn));
        }
    }
}
#endif


TEST(GeneralFrmHam, Elements) {
    const ham_t benchmark = 0.01759459248922075;
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    ASSERT_EQ(frm_ham.m_ints.m_2e->sym(), integrals_2e::syms::DHR);
    Hamiltonian h(&frm_ham);
    {
        buffered::FrmOnv src(h.m_basis);
        src = {{0, 1, 4}, {0, 2, 4}};
        buffered::FrmOnv dst(h.m_basis);
        dst = {{0, 1, 3}, {0, 1, 4}};
        auto helem = h.get_element(src, dst);
        ASSERT_NEAR_EQ(benchmark, helem);
    }
    {
        buffered::FrmBosOnv src(h.m_basis);
        src.m_frm = {{0, 1, 4}, {0, 2, 4}};
        buffered::FrmBosOnv dst(h.m_basis);
        dst.m_frm = {{0, 1, 3}, {0, 1, 4}};
        auto helem = h.get_element(src, dst);
        ASSERT_NEAR_EQ(benchmark, helem);
    }
}

TEST(GeneralFrmHam, RhfEnergy) {
    const ham_comp_t benchmark = -108.76171800006861;
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    Hamiltonian ham(&frm_ham);
    uintv_t chk_orbsyms = {0, 2, 1, 5, 6, 4};
    ASSERT_EQ(ham.m_basis.m_frm.m_abgrp_map.m_site_irreps, chk_orbsyms);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_basis);
    mbf::set_aufbau_mbf(onv, ham.default_particles().m_frm);
    auto elem = ham.get_element(onv);
    ASSERT_NEAR_EQ(elem, benchmark);
    ASSERT_NEAR_EQ(ham.get_energy(onv), benchmark);
}

TEST(GeneralFrmHam, RhfEnergyMolcas) {
    /*
     * HF:      -108.9540866268
     * CASCI:   -109.02180323
     */
    const ham_comp_t benchmark = -108.9540866268;
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5"});
    Hamiltonian ham(&frm_ham);
    uintv_t chk_orbsyms = {0, 0, 0, 0, 0, 0};
    ASSERT_EQ(ham.m_basis.m_frm.m_abgrp_map.m_site_irreps, chk_orbsyms);
    ASSERT_TRUE(ham.m_frm.m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_basis);
    mbf::set_aufbau_mbf(onv, ham.default_particles().m_frm);
    auto elem = ham.get_element(onv);
    ASSERT_NEAR_EQ(elem, benchmark);
    ASSERT_NEAR_EQ(ham.get_energy(onv), benchmark);
}

TEST(GeneralFrmHam, RhfBrillouinTheorem) {
    GeneralFrmHam ham({PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP"});
    ASSERT_TRUE(ham.m_kramers_attrs.conserving());
    buffered::FrmOnv onv(ham.m_basis);
    mbf::set_aufbau_mbf(onv, ham.default_nelec());

    decoded_mbf::frm::SimpleOccs occs(onv);
    decoded_mbf::frm::SimpleVacs vacs(onv);

    buffered::FrmOnv excited(ham.m_basis);
    conn::FrmOnv conn(onv);

    for (uint_t occ : occs.get()){
        for (uint_t vac : vacs.get()){
            conn.m_ann.set(occ);
            conn.m_cre.set(vac);
            ASSERT_NEAR_EQ(ham.get_element_1100(onv, conn), 0.0);
        }
    }
}


TEST(GeneralFrmHam, NonHermitian) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/TC_Be_6-31G/FCIDUMP"});
    ASSERT_EQ(frm_ham.m_ints.m_2e->sym(), integrals_2e::syms::D);
    ASSERT_EQ(frm_ham.m_ints.m_2e->get(5, 2, 1, 6), 0.28747909401816499E-001);
    ASSERT_EQ(frm_ham.m_ints.m_2e->get(2, 5, 6, 1), 0.28747909401816499E-001);
}

TEST(GeneralFrmHam, FromMolcasHdf5Archive) {
    GeneralFrmHam frm_ham({PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5"});
}

TEST(GeneralFrmHam, UhfFcidump) {
    GeneralFrmHam frm_ham_minor({PROJECT_ROOT"/assets/N_UHF/FCIDUMP.spin_minor", FcidumpInfo::SpinMinor});
    GeneralFrmHam frm_ham_major({PROJECT_ROOT"/assets/N_UHF/FCIDUMP.spin_major", FcidumpInfo::SpinMajor});
    GeneralFrmHam frm_ham_blocks({PROJECT_ROOT"/assets/N_UHF/FCIDUMP.spin_blocks", FcidumpInfo::SpinBlocks});
    ASSERT_EQ(frm_ham_minor.m_ints.m_2e->sym(), integrals_2e::syms::DHR);
    ASSERT_EQ(frm_ham_major.m_ints.m_2e->sym(), integrals_2e::syms::DHR);
    ASSERT_EQ(frm_ham_blocks.m_ints.m_2e->sym(), integrals_2e::syms::DHR);

    auto test_fn = [&](uint_t i, uint_t j, uint_t k, uint_t l, ham_t bench_value) {
        ASSERT_NEAR_EQ(frm_ham_minor.m_ints.m_2e->get(i, j, k, l), bench_value);
        ASSERT_NEAR_EQ(frm_ham_major.m_ints.m_2e->get(i, j, k, l), bench_value);
        ASSERT_NEAR_EQ(frm_ham_blocks.m_ints.m_2e->get(i, j, k, l), bench_value);
    };
    /*
     * FCIDUMP.spin_major row
     * 0.03305194106876388    5    2   13   11
     *  < 4 12 | 1 10 >
     */
    test_fn(4, 12, 1, 10,  0.03305194106876388);

    /*
     * FCIDUMP.spin_major row
     * 0.7366168614270663   14   14    1    1
     *  < 13 0 | 13 0 >
     */
    test_fn(13, 0, 13, 0,  0.7366168614270663);
}