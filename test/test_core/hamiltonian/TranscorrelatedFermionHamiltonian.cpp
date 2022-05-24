/**
 * @file TranscorrelatedFermionHamiltonian.cpp
 * @author jph
 * @brief test file for transcorrelated Fermion Hamiltonians
 * @date 2022-05-03
 *
 */

#include <M7_lib/hamiltonian/TcFrmHam.h> // what's being tested
#include <M7_lib/util/consts.h>
#include <gtest/gtest.h>

// #include <iomanip>
// TODO:
// [x] put appropriate TCDUMP(s) into the assets folder
// [ ] get_coeff2200 check for non-Hermiticity (should be handled fine)
// [ ] get_coeff3300 same as get_element3300 up to sign
//          test with and without the parity
// [ ] contracted elements (get_element{00,11,22}00)
//
// [ ] maybe also do a "ui test" as done in TCHINT

#ifdef ENABLE_TCHINT

/**
 * @brief checks if nonhermiticity is handled fine. Constructs a Fermion
 *        Hamiltonian with transcorrelation and checks that it can have
 *        non-Hermitian elements
 *        Note only the 2-body terms are non-Hermitian
 *        [ij|kl] = [ji|lk] from hermiticity
 */
TEST(TranscorrelatedFermionHamiltonian, check_nonhermiticity) {
    // TODO stub
    TcFrmHam ham("FCIDUMP", false, 0);
    // test e.g.
    // [94|12] /= [49|21] (in this case, true)
    auto el1 = ham.get_coeff_2200(9, 4, 1, 2);
    auto el2 = ham.get_coeff_2200(4, 9, 2, 1);
    std::cout << "els:\n" << el1/3.0 << std::endl << el2/3.0 << std::endl;
    ASSERT_FALSE(consts::nearly_equal(el1, el2));
    ASSERT(false);
}

/**
 * @brief checks if get_coeff_element3300 and get_coeff_element3300 are the same
 *        up to parity
 *
 */
TEST(TranscorrelatedFermionHamiltonian, coeff_element3300_parity) {
    // TC Fermion Hamiltonian to be tested
    TcFrmHam ham("FCIDUMP", false, 0);
    buffered::FrmOnv onv(ham.m_nsite);
    onv = {{0, 1},{0, 1}};
    conn::FrmOnv conn(onv);
    // spin-orbital indices to annihilate
    // (one integer -> spin-orbital; pair -> spin, spatial orbital)
    conn.m_cre.add({0, 2});
    conn.m_cre.add({0, 3});
    conn.m_cre.add({1, 9});
    conn.m_ann.add({0, 0}); // alpha electron in spatial orbital 0 to be annihilated
    conn.m_ann.add({0, 1}); // alpha e in spatial orb 1
    conn.m_ann.add({1, 0}); // beta e in spatial orb 0
    auto other = onv;
    conn.apply(onv, other);
    std::cout << other << std::endl;

    defs::ham_t matel33 = ham.get_element_3300(onv, conn);
    // must /3 if the TCDUMP input is in ASCII format
    // expect element (1,2,1|3,4,10)=(1,1,1|3,10,4) (1-based indexing)
    // Also has a phase of -1
    std::cout << "lmat:" << std::endl;
    auto tmp = ham.get_lmat_coeff(1,0,0,3,2,9)/3;
    std::cout << tmp << std::endl;
    //   1 1 2 3 10 4 + 2 1 1 4 3 10 + 1 2 1 10 4 3
    // - 1 1 2 10 3 4 - 1 2 1 4 10 3 - 2 1 1 3 4 10
    // =
    //   1 1 2 3 10 4 + 1 1 2 10 3 4 + 1 1 2 10 3 4
    // - 1 1 2 10 3 4 - 1 1 2 3 4 10 - 1 1 2 4 10 3
    // all but the last two are zero as far as I can tell
    auto benchmark = -0.20360278843803472E-006
                   + 0.0
                   + 0.0
                   - 0.0
                   - 0.0
                   - -0.20360278843803271E-006;
    ASSERT_FLOAT_EQ(-matel33/3.0, benchmark);
    // also test spin-non-conserving excitation (should give 0 for Fermions):
    onv = {{0, 1},{0, 1}};
    // conn::FrmOnv conn(onv);
    // conn.m_ann = {};
    // conn.m_cre = {};
    conn.clear();
    conn.m_ann.add({0, 0});
    conn.m_ann.add({0, 1});
    conn.m_ann.add({1, 0});
    conn.m_cre.add({1, 2});
    conn.m_cre.add({1, 3});
    conn.m_cre.add({1, 9});
    matel33 = ham.get_element_3300(onv, conn);
    // does not conserve spin, ipso facto = 0
    ASSERT_FLOAT_EQ(matel33/3.0, 0.0);

    onv = {{0, 1},{0, 13}};
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
    ASSERT_FLOAT_EQ(matel33/3.0, benchmark);
}

#endif  // ENABLE_TCHINT

