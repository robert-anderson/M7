/**
 * @file TranscorrelatedFermionHamiltonian.cpp
 * @author jph
 * @brief test file for transcorrelated Fermion Hamiltonians
 * @date 2022-05-03
 *
 */

#include <test_core/defs.h>
#include <M7_lib/field/Mbf.h>
#include <M7_lib/hamiltonian/frm/GeneralFrmHam.h>
#include <M7_lib/hamiltonian/frm/TcFrmHam.h>  // what's being tested
#include <M7_lib/util/FpTol.h>
#include <M7_lib/io/Symlink.h>

#ifdef ENABLE_TCHINT

/**
 * @brief check the get_element_0000 method, esp testing contraction
 * (fermion one-particle)
 *
 */
TEST(TranscorrelatedFermionHamiltonian, test_get_element_0000) {
    AssetSymlink tcdump("TC_Be_CCPVDZ/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_CCPVDZ/FCIDUMP", "FCIDUMP");
    // 2.3545134053388529E-002 for 1,2,3,4 (this is only the lmat part)
    // assuming GeneralFrmHam is working properly
    TcFrmHam ham({"FCIDUMP", false});
    GeneralFrmHam two_body_ham({"FCIDUMP", false});
    buffered::FrmOnv onv(ham.m_basis);
    mbf::set_aufbau_mbf(onv, ham.default_nelec());
    auto elem = ham.get_element_0000(onv);
    auto benchmark = two_body_ham.get_element_0000(onv) + 2.3545134053388529E-002;
    ASSERT_NEAR_EQ(elem, benchmark);
}

/**
 * @brief check the get_element_1100 method, esp testing contraction
 * (fermion single excitation)
 *
 */
TEST(TranscorrelatedFermionHamiltonian, test_get_element_1100) {
    AssetSymlink tcdump("TC_Be_CCPVDZ/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_CCPVDZ/FCIDUMP", "FCIDUMP");
    // TCHInt benchmark:
    //      D1=           1           2           3           4
    //  excit1=           4          12
    //   -2.0033789485348489E-003
    TcFrmHam ham({"FCIDUMP", false});
    GeneralFrmHam two_body_ham({"FCIDUMP", false});
    buffered::FrmOnv onv(ham.m_basis);
    mbf::set_aufbau_mbf(onv, ham.default_nelec());
    conn::FrmOnv conn(onv);
    // (one integer -> spin-orbital; pair -> spin, spatial orbital)
    // (spin-minor) spin-orbital 3 goes to spin-orbital 11
    // spinorb 4 -> (orb,spin)=(2,0) -> {0,1}
    conn.m_ann.add({0, 1});
    // spinorb 12 -> (orb,spin)=(6,0) -> {0,5}
    conn.m_cre.add({0, 5});
    auto elemdiff =
        ham.get_element_1100(onv, conn) - two_body_ham.get_element_1100(onv, conn);
    auto benchmark = -2.0033789485348489E-003;
    ASSERT_NEAR_EQ(elemdiff, benchmark);
}

/**
 * @brief check the get_element_2200 method, esp testing contraction
 * (fermion double excitation)
 *
 */
TEST(TranscorrelatedFermionHamiltonian, test_get_element_2200) {
    AssetSymlink tcdump("TC_Be_CCPVDZ/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_CCPVDZ/FCIDUMP", "FCIDUMP");
    // TCHInt benchmark:
    //  D2=           1           2           3           4
    //  excit2=           1           9           4          16
    //   -2.2700965657479885E-005
    // TcFrmHam with spin *minor* ordering, i.e. mirroring NECI & TCHInt
    TcFrmHam ham({"FCIDUMP", false});
    GeneralFrmHam two_body_ham({"FCIDUMP", false});
    buffered::FrmOnv onv(ham.m_basis);
    mbf::set_aufbau_mbf(onv, ham.default_nelec());
    conn::FrmOnv conn(onv);
    // spinorb 4 -> (orb,spin)=(2,0) -> {0, 1}
    conn.m_ann.add({0, 1});
    // spinorb 16 -> (orb,spin)=(8,0) -> {0, 7}
    conn.m_cre.add({0, 7});
    // spinorb 1 -> (orb,spin)=(1,1) -> {1, 0}
    conn.m_ann.add({1, 0});
    // spinorb 9 -> (orb,spin)=(5,1) -> {1, 4}
    conn.m_cre.add({1, 4});
    auto elemdiff =
        ham.get_element_2200(onv, conn) - two_body_ham.get_element_2200(onv, conn);
    auto benchmark = -2.2700965657479885E-005;
    // we also pick up a phase, hence -benchmark
    ASSERT_NEAR_EQ(elemdiff, -benchmark);
}

/**
 * @brief checks if nonhermiticity is handled fine. Constructs a Fermion
 *        Hamiltonian with transcorrelation and checks that it can have
 *        non-Hermitian elements
 *        Note only the 2-body terms are non-Hermitian
 *        [ij|kl] = [ji|lk] from hermiticity
 */
TEST(TranscorrelatedFermionHamiltonian, check_nonhermiticity) {
    AssetSymlink tcdump("TC_Be_CCPVDZ/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_CCPVDZ/FCIDUMP", "FCIDUMP");
    // [ij|kl]=[ji|lk] if Hermitian (chemist notation)
    // remember we antisymmetrise: [ij|kl] - [il|kj]
    // (see FCIDUMP file)
    // the FCIDUMP is in chemist notation but the code is in physicist notations
    TcFrmHam ham({"FCIDUMP", false});
    // these two would be the same assuming Hermiticity, but not in this FCIDUMP
    // chemist notation: 1237 - 1732
    // -0.37788782091129145E-003 - -0.40412978632087910E-003
    auto el1 = ham.get_coeff_2200(0, 2, 1, 6);
    // chemist notation: 2173 - 7123
    // 0.19621434645836822E-002 - -0.63747363394572173E-002
    auto el2 = ham.get_coeff_2200(1, 6, 0, 2);

    ASSERT_NEAR_EQ(el1, el2);
}

/**
 * @brief checks if get_coeff_element3300 and get_coeff_element3300 are the same
 *        up to parity
 *
 */
TEST(TranscorrelatedFermionHamiltonian, coeff_element3300_parity) {
    AssetSymlink tcdump("TC_Be_CCPVDZ/TCDUMP", "TCDUMP");
    AssetSymlink fcidump("TC_Be_CCPVDZ/FCIDUMP", "FCIDUMP");
    // TC Fermion Hamiltonian to be tested
    TcFrmHam ham({"FCIDUMP", false});
    buffered::FrmOnv onv(ham.m_basis);
    onv = {{0, 1}, {0, 1}};
    conn::FrmOnv conn(onv);
    // spin-orbital indices to annihilate
    // (one integer -> spin-orbital; pair -> spin, spatial orbital)
    conn.m_cre.add({0, 2});
    conn.m_cre.add({0, 3});
    conn.m_cre.add({1, 9});
    conn.m_ann.add({0, 0});  // alpha e in spatial orbital 0 to be annihilated
    conn.m_ann.add({0, 1});  // alpha e in spatial orb 1
    conn.m_ann.add({1, 0});  // beta e in spatial orb 0
    auto other = onv;
    conn.apply(onv, other);
    std::cout << other << std::endl;

    ham_t matel33 = ham.get_element_3300(onv, conn);
    // must /3 if the TCDUMP input is in ASCII format
    // expect element (1,2,1|3,4,10)=(1,1,1|3,10,4) (1-based indexing)
    // Also has a phase of -1
    std::cout << "lmat:" << std::endl;
    auto tmp = ham.get_lmat_coeff(1, 0, 0, 3, 2, 9) / 3;
    std::cout << tmp << std::endl;
    // these are the matrix elements that come from anti-symmetrising:
    //   1 1 2 3 10 4 + 2 1 1 4 3 10 + 1 2 1 10 4 3
    // - 1 1 2 10 3 4 - 1 2 1 4 10 3 - 2 1 1 3 4 10
    // =
    //   1 1 2 3 10 4 + 1 1 2 10 3 4 + 1 1 2 10 3 4
    // - 1 1 2 10 3 4 - 1 1 2 3 4 10 - 1 1 2 4 10 3
    // all but the last two are zero
    auto benchmark = -0.20360278843803472E-006 + 0.0 + 0.0 - 0.0 - 0.0 -
                     -0.20360278843803271E-006;
    ASSERT_NEAR_EQ(-matel33 / 3.0, benchmark);
    // also test spin-non-conserving excitation (should give 0 for Fermions):
    onv = {{0, 1}, {0, 1}};
    conn.clear();
    conn.m_ann.add({0, 0});
    conn.m_ann.add({0, 1});
    conn.m_ann.add({1, 0});
    conn.m_cre.add({1, 2});
    conn.m_cre.add({1, 3});
    conn.m_cre.add({1, 9});
    matel33 = ham.get_element_3300(onv, conn);
    // does not conserve spin, ipso facto = 0
    ASSERT_NEAR_EQ(matel33 / 3.0, 0.0);

    onv = {{0, 1}, {0, 13}};
    conn.clear();
    // (one integer -> spin-orbital; pair -> spin, spatial orbital)
    conn.m_cre.add({0, 2});
    conn.m_cre.add({0, 3});
    conn.m_cre.add({1, 9});
    conn.m_ann.add({0, 0});
    conn.m_ann.add({0, 1});
    conn.m_ann.add({1, 0});
    matel33 = ham.get_element_3300(onv, conn);
    // same test as the first one but with positive phase
    ASSERT_NEAR_EQ(matel33 / 3.0, benchmark);
}

#endif  // ENABLE_TCHINT
