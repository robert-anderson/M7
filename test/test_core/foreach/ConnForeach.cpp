//
// Created by anderson on 17/02/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/foreach/ConnForeach.h"


TEST(ConnForeach, FrmGeneralEx1100) {
    const size_t nsite = 8;
    defs::inds setbits = {1, 4, 6, 9, 12};
    buffered::FrmOnv mbf(nsite);
    mbf = setbits;
    auto &clrbits = mbf.m_decoded.m_simple_vacs.get();

    size_t iiter = 0ul;
    auto fn = [&](const conn::FrmOnv &conn) {
        const auto cre = conn.m_cre[0];
        const auto ann = conn.m_ann[0];
        const auto iann = iiter / clrbits.size();
        const auto icre = iiter - iann * clrbits.size();
        ASSERT_EQ(cre, clrbits[icre]);
        ASSERT_EQ(ann, setbits[iann]);
        ++iiter;
    };
    conn_foreach::frm::General<1> foreach(nsite);
    ASSERT_EQ(foreach.m_exsig, exsig_utils::ex_single);
    foreach.loop_fn(mbf, fn);
    ASSERT_EQ(iiter, setbits.size() * clrbits.size());
}

TEST(ConnForeach, FrmGeneralEx2200) {
    const size_t nsite = 4;
    defs::inds setbits = {1, 3, 5, 6};
    buffered::FrmOnv mbf(nsite);
    mbf = setbits;
    defs::inds clrbits = {0, 2, 4, 7};
    ASSERT_EQ(mbf.m_decoded.m_simple_vacs.get(), clrbits);
    std::vector<defs::inds> setbit_pairs = {
            {1, 3},
            {1, 5},
            {3, 5},
            {1, 6},
            {3, 6},
            {5, 6}};
    std::vector<defs::inds> clrbit_pairs = {
            {0, 2},
            {0, 4},
            {2, 4},
            {0, 7},
            {2, 7},
            {4, 7}};

    size_t iiter = 0ul;
    auto fn = [&](const conn::FrmOnv &conn) {
        const auto iann_pair = iiter / clrbit_pairs.size();
        const auto icre_pair = iiter - iann_pair * clrbit_pairs.size();
        ASSERT_EQ(conn.m_cre, clrbit_pairs[icre_pair]);
        ASSERT_EQ(conn.m_ann, setbit_pairs[iann_pair]);
        ++iiter;
    };
    conn_foreach::frm::General<2> foreach(nsite);
    ASSERT_EQ(foreach.m_exsig, exsig_utils::ex_double);
    foreach.loop_fn(mbf, fn);
    ASSERT_EQ(iiter, setbit_pairs.size() * clrbit_pairs.size());
}

TEST(ConnForeach, FrmMs2ConserveEx1100) {
    const size_t nsite = 8;
    defs::inds alpha_setbits = {1, 4, 6};
    defs::inds beta_setbits = {3, 4, 7};
    buffered::FrmOnv mbf(nsite);
    mbf = {alpha_setbits, beta_setbits};
    const auto &occs = mbf.m_decoded.m_spin_occs.get();
    const auto &vacs = mbf.m_decoded.m_spin_vacs.get();

    const auto nalpha_setbit = alpha_setbits.size();
    const auto nalpha_clrbit = nsite - nalpha_setbit;
    const auto nalpha_pair = nalpha_setbit * nalpha_clrbit;

    const auto nbeta_setbit = beta_setbits.size();
    const auto nbeta_clrbit = nsite - nbeta_setbit;
    const auto nbeta_pair = nbeta_setbit * nbeta_clrbit;

    size_t iiter = 0ul;
    auto fn = [&](const conn::FrmOnv &conn) {
        if (iiter < nalpha_pair) {
            // a -> a sector
            const auto iann = iiter / nalpha_clrbit;
            const auto icre = iiter - iann * nalpha_clrbit;
            ASSERT_EQ(occs[0][iann], conn.m_ann[0]);
            ASSERT_EQ(vacs[0][icre], conn.m_cre[0]);
        } else {
            const auto iiter_beta = iiter - nalpha_pair;
            // b -> b sector
            const auto iann = iiter_beta / nbeta_clrbit;
            const auto icre = iiter_beta - iann * nbeta_clrbit;
            ASSERT_EQ(occs[1][iann], conn.m_ann[0]);
            ASSERT_EQ(vacs[1][icre], conn.m_cre[0]);
        }
        ++iiter;
    };
    conn_foreach::frm::Ms2Conserve<1> foreach(nsite);
    ASSERT_EQ(foreach.m_exsig, exsig_utils::ex_single);
    foreach.loop_fn(mbf, fn);
    ASSERT_EQ(iiter, nalpha_pair + nbeta_pair);
}


TEST(ConnForeach, FrmMs2ConserveEx2200) {
    const size_t nsite = 8;
    defs::inds alpha_setbits = {1, 2, 5};
    defs::inds beta_setbits = {0, 2, 4, 6};
    buffered::FrmOnv mbf(nsite);
    mbf = {alpha_setbits, beta_setbits};

    const auto nalpha_setbit = alpha_setbits.size();
    const auto npair_alpha_setbit = integer_utils::nspair(nalpha_setbit);
    const auto nalpha_clrbit = nsite - nalpha_setbit;
    const auto npair_alpha_clrbit = integer_utils::nspair(nalpha_clrbit);

    const auto nbeta_setbit = beta_setbits.size();
    const auto npair_beta_setbit = integer_utils::nspair(nbeta_setbit);
    const auto nbeta_clrbit = nsite - nbeta_setbit;
    const auto npair_beta_clrbit = integer_utils::nspair(nbeta_clrbit);

    size_t iiter = 0ul;

    const size_t naaaa = npair_alpha_clrbit * npair_alpha_setbit;
    const size_t nabab = nalpha_setbit * nalpha_clrbit * nbeta_setbit * nbeta_clrbit;
    const size_t nbbbb = npair_beta_clrbit * npair_beta_setbit;

    auto fn = [&](const conn::FrmOnv &conn) {
        if (iiter < naaaa) {
            // aa -> aa
            ASSERT_EQ(conn.m_cre.nalpha(), 2);
            ASSERT_EQ(conn.m_ann.nalpha(), 2);
        } else if (iiter < naaaa + nabab) {
            // ab -> ab
            ASSERT_EQ(conn.m_cre.nalpha(), 1);
            ASSERT_EQ(conn.m_ann.nalpha(), 1);
        } else {
            // bb -> bb
            ASSERT_EQ(conn.m_cre.nalpha(), 0);
            ASSERT_EQ(conn.m_ann.nalpha(), 0);
        }
        ++iiter;
    };
    conn_foreach::frm::Ms2Conserve<2> foreach(nsite);
    ASSERT_EQ(foreach.m_exsig, exsig_utils::ex_double);
    foreach.loop_fn(mbf, fn);
    ASSERT_EQ(iiter, naaaa + nabab + nbbbb);
}
