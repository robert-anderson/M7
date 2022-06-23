//
// Created by Robert J. Anderson on 17/02/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/foreach/ConnForeach.h"

namespace conn_foreach_test {
    struct Result {
        const defs::ivec_t m_ann, m_cre;

        Result(defs::ivec_t ann, defs::ivec_t cre) : m_ann(std::move(ann)), m_cre(std::move(cre)) {}

        Result(size_t ann, size_t cre) : m_ann({ann}), m_cre({cre}) {}
    };

    typedef std::vector<Result> results_t;

    results_t product_results(const std::vector<defs::ivec_t> &anns, const std::vector<defs::ivec_t> &cres) {
        const std::vector<defs::ivec_t> default_ops = {{}};
        const auto &anns_ref = anns.empty() ? default_ops : anns;
        const auto &cres_ref = cres.empty() ? default_ops : cres;

        results_t results;
        results.reserve(anns_ref.size() * cres_ref.size());
        for (auto &ann: anns_ref) {
            for (auto &cre: cres_ref) {
                results.push_back({ann, cre});
            }
        }
        return results;
    }

    static defs::ivec_t creatable_mode_indices(const field::BosOnv &mbf, size_t nboson_max) {
        defs::ivec_t inds;
        for (size_t imode = 0ul; imode < mbf.m_basis.m_nmode; ++imode) {
            size_t nocc = mbf[imode];
            if (nocc < nboson_max) inds.push_back(imode);
        }
        return inds;
    }

    static defs::ivec_t annihilatable_mode_indices(const field::BosOnv &mbf) {
        defs::ivec_t inds;
        for (size_t imode = 0ul; imode < mbf.m_basis.m_nmode; ++imode) if (mbf[imode]) inds.push_back(imode);
        return inds;
    }
}

TEST(ConnForeach, FrmGeneralEx1100FrmOnv) {
    const size_t nsite = 8;
    defs::ivec_t setbits = {1, 4, 6, 9, 12};
    buffered::FrmOnv mbf(nsite);
    mbf = setbits;
    auto &clrbits = mbf.m_decoded.m_simple_vacs.get();

    conn::FrmOnv conn(mbf.m_basis);
    size_t iiter = 0ul;
    auto fn = [&]() {
        const auto cre = conn.m_cre[0];
        const auto ann = conn.m_ann[0];
        const auto iann = iiter / clrbits.size();
        const auto icre = iiter - iann * clrbits.size();
        ASSERT_EQ(cre, clrbits[icre]);
        ASSERT_EQ(ann, setbits[iann]);
        ++iiter;
    };
    conn_foreach::frm::General<1> foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_single);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(iiter, setbits.size() * clrbits.size());
    // then, the run time polymorphic loop:
    iiter = 0;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, setbits.size() * clrbits.size());

}

TEST(ConnForeach, FrmGeneralEx1100FrmBosOnv) {
    const size_t nsite = 8;
    defs::ivec_t setbits = {1, 4, 6, 9, 12};
    buffered::FrmBosOnv mbf(nsite, 0ul);
    mbf.m_frm = setbits;
    auto &clrbits = mbf.m_frm.m_decoded.m_simple_vacs.get();
    conn::FrmBosOnv conn(mbf);

    size_t iiter = 0ul;
    auto fn = [&]() {
        const auto cre = conn.m_frm.m_cre[0];
        const auto ann = conn.m_frm.m_ann[0];
        const auto iann = iiter / clrbits.size();
        const auto icre = iiter - iann * clrbits.size();
        ASSERT_EQ(cre, clrbits[icre]);
        ASSERT_EQ(ann, setbits[iann]);
        ++iiter;
    };
    conn_foreach::frm::General<1> foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_single);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn.m_frm, mbf.m_frm, fn);
    ASSERT_EQ(iiter, setbits.size() * clrbits.size());
    // then, the run time polymorphic loop:
    iiter = 0;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, setbits.size() * clrbits.size());
}

TEST(ConnForeach, FrmGeneralEx2200) {
    const size_t nsite = 4;
    defs::ivec_t setbits = {1, 3, 5, 6};
    buffered::FrmOnv mbf(nsite);
    mbf = setbits;
    defs::ivec_t clrbits = {0, 2, 4, 7};
    ASSERT_EQ(mbf.m_decoded.m_simple_vacs.get(), clrbits);
    std::vector<defs::ivec_t> setbit_pairs = {{1, 3},
                                              {1, 5},
                                              {3, 5},
                                              {1, 6},
                                              {3, 6},
                                              {5, 6}};
    std::vector<defs::ivec_t> clrbit_pairs = {{0, 2},
                                              {0, 4},
                                              {2, 4},
                                              {0, 7},
                                              {2, 7},
                                              {4, 7}};
    auto results = conn_foreach_test::product_results(setbit_pairs, clrbit_pairs);

    conn::FrmOnv conn(mbf);
    auto result = results.cbegin();
    auto fn = [&]() {
        ASSERT_EQ(conn.m_ann, result->m_ann);
        ASSERT_EQ(conn.m_cre, result->m_cre);
        ++result;
    };
    conn_foreach::frm::General<2> foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_double);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(result, results.cend());
    // then, the run time polymorphic loop:
    result = results.cbegin();
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(result, results.cend());
}

TEST(ConnForeach, FrmHubbardEx1100) {
    /*
     * 2d hubbard model example det      site indices (irow*4+icol)
     *  x o x o                          0 1 2 3
     *  o o o x                          4 5 6 7
     *  x x o o                          8 ...
     *  periodic BCs top to bottom (major dimension: "row")
     *  open BCs left to right (minor dimension: "col")
     */
    const size_t nrow = 3, ncol = 4;
    const sys::frm::Basis basis(lattice::make("ortho", {nrow, ncol}, {1, 0}));
    // horizontal neighbors
    ASSERT_TRUE(basis.m_lattice->phase(1, 2));
    ASSERT_TRUE(basis.m_lattice->phase(2, 1));
    // vertical neighbors
    ASSERT_TRUE(basis.m_lattice->phase(1, 5));
    ASSERT_TRUE(basis.m_lattice->phase(5, 1));
    // periodic vertical neighbors
    ASSERT_TRUE(basis.m_lattice->phase(1, 9));
    ASSERT_TRUE(basis.m_lattice->phase(9, 1));
    // horizontal neighbors
    ASSERT_TRUE(basis.m_lattice->phase(4, 5));
    ASSERT_TRUE(basis.m_lattice->phase(5, 4));
    // NOT periodic horizontal neighbors
    ASSERT_FALSE(basis.m_lattice->phase(4, 7));
    ASSERT_FALSE(basis.m_lattice->phase(7, 4));

    buffered::FrmOnv mbf(basis);
    mbf = {0, 2, 7, 8, 9};

    conn_foreach_test::results_t results = {
            {0, 4},
            {0, 1},
            {2, 10},
            {2, 6},
            {2, 1},
            {2, 3},
            {7, 3},
            {7, 11},
            {7, 6},
            {8, 4},
            {9, 5},
            {9, 1},
            {9, 10}};
    conn::FrmOnv conn(basis);
    auto result = results.cbegin();
    auto fn = [&]() {
        ASSERT_EQ(conn.m_cre, result->m_cre);
        ASSERT_EQ(conn.m_ann, result->m_ann);
        ++result;
    };
    conn_foreach::frm::Hubbard foreach;
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(result, results.cend());
    // then, the run time polymorphic loop:
    result = results.cbegin();
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(result, results.cend());
}

TEST(ConnForeach, FrmHeisenbergEx2200) {
    /*
     * 2d heisenberg model example alphas
     *  1 0 0 1 0 1 1 0
     *  periodic BCs
     */
    const size_t nsite = 8;
    const sys::frm::Basis basis(lattice::make("ortho", {nsite}, {1}));
    buffered::FrmOnv mbf(basis);
    mbf.set_spins({0, 3, 5});

    conn_foreach_test::results_t results = {
            {{0, 15}, {7, 8}},
            {{0, 9},  {1, 8}},
            {{3, 10}, {2, 11}},
            {{3, 12}, {4, 11}},
            {{5, 12}, {4, 13}},
            {{5, 14}, {6, 13}}};

    auto result = results.cbegin();
    conn::FrmOnv conn(mbf);
    auto fn = [&]() {
        ASSERT_EQ(conn.exsig(), exsig::ex_double);
        ASSERT_EQ(conn.m_cre, result->m_cre);
        ASSERT_EQ(conn.m_ann, result->m_ann);
        ++result;
    };
    conn_foreach::frm::Heisenberg foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_double);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(result, results.cend());
    // then, the run time polymorphic loop:
    result = results.cbegin();
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(result, results.cend());
}

TEST(ConnForeach, FrmMs2ConserveEx1100) {
    const size_t nsite = 8;
    defs::ivec_t alpha_setbits = {1, 4, 6};
    defs::ivec_t beta_setbits = {3, 4, 7};
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
    conn::FrmOnv conn(mbf);
    auto fn = [&]() {
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
    conn_foreach::frm::Ms2Conserve<1> foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_single);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(iiter, nalpha_pair + nbeta_pair);
    // then, the run time polymorphic loop:
    iiter = 0;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, nalpha_pair + nbeta_pair);
}

TEST(ConnForeach, FrmMs2ConserveEx2200) {
    const size_t nsite = 8;
    defs::ivec_t alpha_setbits = {1, 2, 5};
    defs::ivec_t beta_setbits = {0, 2, 4, 6};
    buffered::FrmOnv mbf(nsite);
    mbf = {alpha_setbits, beta_setbits};

    const auto nalpha_setbit = alpha_setbits.size();
    const auto npair_alpha_setbit = integer::nspair(nalpha_setbit);
    const auto nalpha_clrbit = nsite - nalpha_setbit;
    const auto npair_alpha_clrbit = integer::nspair(nalpha_clrbit);

    const auto nbeta_setbit = beta_setbits.size();
    const auto npair_beta_setbit = integer::nspair(nbeta_setbit);
    const auto nbeta_clrbit = nsite - nbeta_setbit;
    const auto npair_beta_clrbit = integer::nspair(nbeta_clrbit);

    size_t iiter = 0ul;

    const size_t naaaa = npair_alpha_clrbit * npair_alpha_setbit;
    const size_t nabab = nalpha_setbit * nalpha_clrbit * nbeta_setbit * nbeta_clrbit;
    const size_t nbbbb = npair_beta_clrbit * npair_beta_setbit;

    conn::FrmOnv conn(mbf);
    auto fn = [&]() {
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
    conn_foreach::frm::Ms2Conserve<2> foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_double);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(iiter, naaaa + nabab + nbbbb);
    // then, the run time polymorphic loop:
    iiter = 0ul;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, naaaa + nabab + nbbbb);
}

TEST(ConnForeach, BosEx0001) {
    const size_t nmode = 6;
    defs::ivec_t occs = {0, 2, 0, 1, 5, 1};
    buffered::BosOnv mbf(nmode);
    mbf = occs;
    auto &chk_modes = mbf.m_decoded.m_occ_modes.get();

    auto iiter = 0ul;
    conn::BosOnv conn(mbf);
    auto fn = [&]() {
        ASSERT_EQ(conn.m_ann[0].m_imode, chk_modes[iiter]);
        ASSERT_FALSE(conn.m_cre.size());
        ++iiter;
    };
    conn_foreach::bos::Ann foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_0001);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(iiter, chk_modes.size());
    // then, the run time polymorphic loop:
    iiter = 0ul;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, chk_modes.size());
}

TEST(ConnForeach, BosEx0010BosOnv) {
    const size_t nmode = 6;
    defs::ivec_t occs = {0, 2, 0, 1, 5, 1};

    {
        const size_t occ_cutoff = 10;
        buffered::BosOnv mbf(nmode, occ_cutoff);
        mbf = occs;
        defs::ivec_t chk_modes = {0, 1, 2, 3, 4, 5};
        auto iiter = 0ul;
        conn::BosOnv conn(mbf);
        auto fn = [&]() {
            ASSERT_FALSE(conn.m_ann.size());
            ASSERT_EQ(conn.m_cre[0].m_imode, chk_modes[iiter]);
            ++iiter;
        };
        conn_foreach::bos::Cre foreach;
        ASSERT_EQ(foreach.m_exsig, exsig::ex_0010);
        // first, the compile time polymorphic loop:
        foreach.loop_fn(conn, mbf, fn);
        ASSERT_EQ(iiter, chk_modes.size());
        // then, the run time polymorphic loop:
        iiter = 0ul;
        foreach.loop(conn, mbf, fn);
        ASSERT_EQ(iiter, chk_modes.size());
    }

    {
        const size_t occ_cutoff = 5;
        // lower maximum occupation
        buffered::BosOnv mbf(nmode, occ_cutoff);
        mbf = occs;
        defs::ivec_t chk_modes = {0, 1, 2, 3, 5};
        auto iiter = 0ul;
        conn::BosOnv conn(mbf);
        auto fn = [&]() {
            ASSERT_FALSE(conn.m_ann.size());
            ASSERT_EQ(conn.m_cre[0].m_imode, chk_modes[iiter]);
            ++iiter;
        };
        conn_foreach::bos::Cre foreach;
        ASSERT_EQ(foreach.m_exsig, exsig::ex_0010);
        // first, the compile time polymorphic loop:
        foreach.loop_fn(conn, mbf, fn);
        ASSERT_EQ(iiter, chk_modes.size());
        // then, the run time polymorphic loop:
        iiter = 0ul;
        foreach.loop(conn, mbf, fn);
        ASSERT_EQ(iiter, chk_modes.size());
    }
}

TEST(ConnForeach, BosEx0010FrmBosOnv) {
    const size_t nmode = 6;
    const size_t bos_occ_cutoff = 5;
    defs::ivec_t occs = {0, 2, 0, 1, 5, 1};
    // lower maximum occupation
    buffered::FrmBosOnv mbf(0ul, nmode, bos_occ_cutoff);
    mbf.m_bos = occs;
    defs::ivec_t chk_modes = {0, 1, 2, 3, 5};
    auto iiter = 0ul;
    conn::FrmBosOnv conn(mbf);
    auto fn = [&]() {
        ASSERT_FALSE(conn.m_bos.m_ann.size());
        ASSERT_EQ(conn.m_bos.m_cre[0].m_imode, chk_modes[iiter]);
        ++iiter;
    };
    conn_foreach::bos::Cre foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_0010);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn.m_bos, mbf.m_bos, fn);
    ASSERT_EQ(iiter, chk_modes.size());
    // then, the run time polymorphic loop:
    iiter = 0ul;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, chk_modes.size());
}

TEST(ConnForeach, FrmBosEx1110) {
    const size_t nsite = 6, nmode = 8;
    for (size_t bos_occ_cutoff = 0ul; bos_occ_cutoff < 6; ++bos_occ_cutoff) {
        buffered::FrmBosOnv mbf(nsite, nmode, bos_occ_cutoff);

        ASSERT_EQ(mbf.m_bos.m_basis.m_nmode, nmode);
        mbf = {{0, 4, 6, 8, 11},
               {2, 0, 1, 0, 1, 4, 0, 5}};
        using namespace conn_foreach;
        /*
         * first, do the product manually. The two components of this product iterator are already tested above, here we
         * only need to ensure that the nested loop is performed as required
         */
        conn_foreach_test::results_t frm_results;
        {
            conn::FrmOnv conn(mbf.m_frm);
            auto fn = [&]() {
                frm_results.emplace_back(conn.m_ann.inds(), conn.m_cre.inds());
            };
            frm::Ms2Conserve<1> foreach;
            foreach.loop_fn(conn, mbf.m_frm, fn);
        }

        auto mode_inds = conn_foreach_test::creatable_mode_indices(mbf.m_bos, bos_occ_cutoff);
        size_t iiter = 0ul;
        conn::FrmBosOnv conn(mbf);
        auto fn = [&]() {
            auto ifrm_result = iiter / mode_inds.size();
            ASSERT_EQ(conn.m_frm.m_ann, frm_results[ifrm_result].m_ann);
            ASSERT_EQ(conn.m_frm.m_cre, frm_results[ifrm_result].m_cre);
            auto imode = iiter - ifrm_result * mode_inds.size();
            ASSERT_EQ(conn.m_bos.m_cre[0].m_imode, mode_inds[imode]);
            ++iiter;
        };
        frmbos::Product<frm::Ms2Conserve<1>, bos::Cre> foreach;
        ASSERT_EQ(foreach.m_exsig, exsig::ex_1110);
        // first, the compile time polymorphic loop:
        foreach.loop_fn(conn, mbf, fn);
        ASSERT_EQ(iiter, frm_results.size() * mode_inds.size());
        // then, the run time polymorphic loop:
        iiter = 0ul;
        foreach.loop(conn, mbf, fn);
        ASSERT_EQ(iiter, frm_results.size() * mode_inds.size());
    }
}

TEST(ConnForeach, FrmBosEx1101) {
    const size_t nsite = 6, nmode = 8;
    buffered::FrmBosOnv mbf(nsite, nmode);
    mbf = {{0, 4, 6, 8, 11},
           {2, 0, 1, 0, 1, 4, 0, 5}};
    using namespace conn_foreach;
    /*
     * first, do the product manually. The two components of this product iterator are already tested above, here we
     * only need to ensure that the nested loop is performed as required
     */
    conn_foreach_test::results_t frm_results;
    {
        conn::FrmOnv conn(mbf.m_frm);
        auto fn = [&]() {
            frm_results.emplace_back(conn.m_ann.inds(), conn.m_cre.inds());
        };
        frm::Ms2Conserve<1> foreach;
        foreach.loop_fn(conn, mbf.m_frm, fn);
    }

    auto mode_inds = conn_foreach_test::annihilatable_mode_indices(mbf.m_bos);
    size_t iiter = 0ul;
    conn::FrmBosOnv conn(mbf);
    auto fn = [&]() {
        auto ifrm_result = iiter / mode_inds.size();
        ASSERT_EQ(conn.m_frm.m_ann, frm_results[ifrm_result].m_ann);
        ASSERT_EQ(conn.m_frm.m_cre, frm_results[ifrm_result].m_cre);
        auto imode = iiter - ifrm_result * mode_inds.size();
        ASSERT_EQ(conn.m_bos.m_ann[0].m_imode, mode_inds[imode]);
        ++iiter;
    };
    frmbos::Product<frm::Ms2Conserve<1>, bos::Ann> foreach;
    ASSERT_EQ(foreach.m_exsig, exsig::ex_1101);
    // first, the compile time polymorphic loop:
    foreach.loop_fn(conn, mbf, fn);
    ASSERT_EQ(iiter, frm_results.size() * mode_inds.size());
    // then, the run time polymorphic loop:
    iiter = 0ul;
    foreach.loop(conn, mbf, fn);
    ASSERT_EQ(iiter, frm_results.size() * mode_inds.size());
}