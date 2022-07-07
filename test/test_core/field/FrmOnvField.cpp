//
// Created by Robert J. Anderson on 09/02/2022.
//

#include <numeric>

#include "gtest/gtest.h"
#include <M7_lib/foreach/BasicForeach.h>
#include "M7_lib/table/BufferedFields.h"


namespace frm_onv_field_test {
    typedef const std::function<bool(bool, bool)>& include_fn_t;
    typedef const std::function<uint_t(const field::FrmOnv&)>& count_fn_t;
    template<typename foreach_fn_t>
    void site_foreach(include_fn_t include_fn, count_fn_t count_fn, foreach_fn_t foreach_fn,
                      bool fill_alpha= false, bool fill_beta= false) {
        auto test_fn = [&](uint_t nsite, uint_t nset) {
            buffered::FrmOnv mbf(nsite);
            uint_t nsetbit_chk = 0ul;
            uintv_t alpha_setbits(nsite);
            if (fill_alpha) {
                std::iota(alpha_setbits.begin(), alpha_setbits.end(), 0);
                nsetbit_chk+=nsite;
            }
            else {
                alpha_setbits = hash::unique_in_range<uint_t>(0, nset, 0, nsite, true);
                nsetbit_chk+=nset;
            }
            uintv_t beta_setbits(nsite);
            if (fill_beta) {
                std::iota(beta_setbits.begin(), beta_setbits.end(), 0);
                nsetbit_chk+=nsite;
            }
            else {
                beta_setbits = hash::unique_in_range<uint_t>(1, nset, 0, nsite, true);
                nsetbit_chk+=nset;
            }

            mbf = {alpha_setbits, beta_setbits};
            ASSERT_EQ(mbf.nsetbit(), nsetbit_chk);
            uintv_t isites_included;
            for (uint_t isite=0ul; isite<nsite; ++isite){
                if (include_fn(mbf.get({0, isite}), mbf.get({1, isite})))
                    isites_included.push_back(isite);
            }
            ASSERT_EQ(isites_included.size(), count_fn(mbf));
            /*
             * check that the iterator returns the correct site indices
             */
            auto it = isites_included.cbegin();
            const auto it_end = isites_included.cend();
            auto fn = [&it, &it_end](uint_t isite) {
                ASSERT_NE(it, it_end);
                ASSERT_EQ(isite, *(it++));
            };
            foreach_fn(mbf, fn);
            ASSERT_EQ(it, isites_included.cend());
        };
        // less than one dataword total
        test_fn(24, 16);
        // less than one dataword per channel, unaligned
        test_fn(39, 19);
        // two datawords per channel, unaligned
        test_fn(68, 17);
        // multiple datawords per channel, unaligned
        test_fn(235, 145);
        // multiple datawords per channel, aligned
        test_fn(3 * Buffer::c_nbit_word, Buffer::c_nbit_word);
    }
}


TEST(FrmOnvField, SetFromInds) {
    const uint_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    uintv_t setbits = {1, 90, nsite - 1, nsite, 2 * nsite - 1};
    mbf = setbits;
    for (uint_t ibit=0ul; ibit<mbf.m_basis.m_nspinorb; ++ibit){
        bool is_set = std::find(setbits.cbegin(), setbits.cend(), ibit)!=setbits.cend();
        ASSERT_EQ(mbf.get(ibit), is_set);
    }
}

TEST(FrmOnvField, ClrSpinChannel) {
    const uint_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    ASSERT_EQ(mbf.nalpha(), 0ul);
    ASSERT_EQ(mbf.nsetbit(), 0ul);
    mbf.put_spin_channel(0, true);
    ASSERT_EQ(mbf.nalpha(), nsite);
    ASSERT_EQ(mbf.nsetbit(), nsite);
    mbf.put_spin_channel(1, true);
    ASSERT_EQ(mbf.nalpha(), nsite);
    ASSERT_EQ(mbf.nsetbit(), 2*nsite);
    mbf.put_spin_channel(0, false);
    ASSERT_EQ(mbf.nalpha(), 0ul);
    ASSERT_EQ(mbf.nsetbit(), nsite);
    mbf.put_spin_channel(1, false);
    ASSERT_EQ(mbf.nalpha(), 0ul);
    ASSERT_EQ(mbf.nsetbit(), 0ul);
}

TEST(FrmOnvField, ForeachSetBit) {
    const uint_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hash::unique_in_range<uint_t>(0, 64, 0, mbf.m_basis.m_nspinorb, true);
    mbf = setbits;

    auto it = setbits.cbegin();
    auto fn = [&it](uint_t ibit){
        ASSERT_EQ(ibit, *(it++));
    };
    mbf.foreach_setbit(fn);
    // make sure all bits were iterated over
    ASSERT_EQ(it, setbits.cend());
}

TEST(FrmOnvField, ForeachSetBitPair) {
    const uint_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hash::unique_in_range<uint_t>(0, 64, 0, mbf.m_basis.m_nspinorb, true);
    mbf = setbits;

    /*
     * set up checking vectors for the pairs of set bits that should be found
     * the order of set bit pairs returned should match that of the ascending-ordered basic foreach pair iterator
     */
    using namespace basic_foreach;
    v_t<ctnd::inds_t<2>> setbit_pairs;
    {
        auto fn = [&](const ctnd::inds_t<2>& inds) {
            setbit_pairs.push_back({setbits[inds[0]], setbits[inds[1]]});
        };
        ctnd::Ordered<2, true, true> foreach(setbits.size());
        foreach.loop(fn);
    }

    auto it = setbit_pairs.cbegin();
    auto fn = [&it](uint_t ibit, uint_t jbit){
        ASSERT_EQ(ibit, (*it)[0]);
        ASSERT_EQ(jbit, (*it)[1]);
        ++it;
    };
    mbf.foreach_setbit_pair(fn);
    // make sure all bits were iterated over
    ASSERT_EQ(it, setbit_pairs.cend());
}

TEST(FrmOnvField, ForeachSetBitTriple) {
    const uint_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hash::unique_in_range<uint_t>(0, 64, 0, mbf.m_basis.m_nspinorb, true);
    mbf = setbits;

    /*
     * set up checking vectors for the pairs of set bits that should be found
     * the order of set bit pairs returned should match that of the ascending-ordered basic foreach pair iterator
     */
    using namespace basic_foreach;
    v_t<ctnd::inds_t<3>> setbit_triples;
    {
        auto fn = [&](const ctnd::inds_t<3>& inds) {
            setbit_triples.push_back({setbits[inds[0]], setbits[inds[1]], setbits[inds[2]]});
        };
        ctnd::Ordered<3, true, true> foreach(setbits.size());
        foreach.loop(fn);
    }

    auto it = setbit_triples.cbegin();
    auto fn = [&it](uint_t ibit, uint_t jbit, uint_t kbit){
        ASSERT_EQ(ibit, (*it)[0]);
        ASSERT_EQ(jbit, (*it)[1]);
        ASSERT_EQ(kbit, (*it)[2]);
        ++it;
    };
    mbf.foreach_setbit_triple(fn);
    // make sure all bits were iterated over
    ASSERT_EQ(it, setbit_triples.cend());
}


TEST(FrmOnvField, ForeachAlphaSetBit) {
    auto include_fn = [](bool alpha, bool /*beta*/){
        return alpha;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nalpha();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_alpha(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn, false, true);
}

TEST(FrmOnvField, ForeachBetaSetBit) {
    auto include_fn = [](bool /*alpha*/, bool beta){
        return beta;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nbeta();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_beta(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn, true, false);
}


TEST(FrmOnvField, ForeachOpenShell) {
    auto include_fn = [](bool alpha, bool beta){
        return alpha != beta;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nopen_shell();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_open_shell(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn);
}

TEST(FrmOnvField, ForeachOpenShellAlpha) {
    auto include_fn = [](bool alpha, bool beta){
        return alpha && !beta;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nopen_shell_alpha();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_open_shell_alpha(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn);
}

TEST(FrmOnvField, ForeachOpenShellBeta) {
    auto include_fn = [](bool alpha, bool beta){
        return !alpha && beta;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nopen_shell_beta();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_open_shell_beta(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn);
}

TEST(FrmOnvField, ForeachClosedShell) {
    auto include_fn = [](bool alpha, bool beta){
        return alpha && beta;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nclosed_shell();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_closed_shell(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn);
}

TEST(FrmOnvField, ForeachOccupiedSite) {
    auto include_fn = [](bool alpha, bool beta){
        return alpha || beta;
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.noccupied_site();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_occupied_site(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn);
}

TEST(FrmOnvField, ForeachUnoccupiedSite) {
    auto include_fn = [](bool alpha, bool beta){
        return !(alpha || beta);
    };
    auto count_fn = [](const field::FrmOnv& mbf){
        return mbf.nunoccupied_site();
    };
    auto test_fn = [](const field::FrmOnv& mbf, std::function<void(uint_t)> body_fn) {
        mbf.foreach_unoccupied_site(body_fn);
    };
    frm_onv_field_test::site_foreach(include_fn, count_fn, test_fn);
}