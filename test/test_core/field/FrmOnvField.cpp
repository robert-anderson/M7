//
// Created by anderson on 09/02/2022.
//

#include "gtest/gtest.h"

#include <M7_lib/foreach/BasicForeach.h>
#include "M7_lib/table/BufferedFields.h"

TEST(FrmOnvField, SetFromInds) {
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    defs::inds setbits = {1, 90, nsite-1, nsite, 2*nsite-1};
    mbf = setbits;
    for (size_t ibit=0ul; ibit<mbf.m_nspinorb; ++ibit){
        bool is_set = std::find(setbits.cbegin(), setbits.cend(), ibit)!=setbits.cend();
        ASSERT_EQ(mbf.get(ibit), is_set);
    }
}

TEST(FrmOnvField, ClrSpinChannel) {
    const size_t nsite = 123;
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
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hashing::unique_in_range(0, 64, 0, mbf.m_nspinorb, true);
    mbf = setbits;

    auto it = setbits.cbegin();
    auto fn = [&it](size_t ibit){
        ASSERT_EQ(ibit, *(it++));
    };
    mbf.foreach_setbit(fn);
    // make sure all bits were iterated over
    ASSERT_EQ(it, setbits.cend());
}

TEST(FrmOnvField, ForeachSetBitPair) {
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hashing::unique_in_range(0, 64, 0, mbf.m_nspinorb, true);
    mbf = setbits;

    /*
     * set up checking vectors for the pairs of set bits that should be found
     * the order of set bit pairs returned should match that of the ascending-ordered basic foreach pair iterator
     */
    using namespace basic_foreach;
    std::vector<ctnd::inds_t<2>> setbit_pairs;
    {
        auto fn = [&](const ctnd::inds_t<2>& inds) {
            setbit_pairs.push_back({setbits[inds[0]], setbits[inds[1]]});
        };
        ctnd::Ordered<2, true, true> foreach(setbits.size());
        foreach.loop(fn);
    }

    auto it = setbit_pairs.cbegin();
    auto fn = [&it](size_t ibit, size_t jbit){
        ASSERT_EQ(ibit, (*it)[0]);
        ASSERT_EQ(jbit, (*it)[1]);
        ++it;
    };
    mbf.foreach_setbit_pair(fn);
    // make sure all bits were iterated over
    ASSERT_EQ(it, setbit_pairs.cend());
}


TEST(FrmOnvField, ForeachOpenShell) {
    const size_t nsite = 68; //123;
    const size_t nset = 17; //64;
    buffered::FrmOnv mbf(nsite);
    auto alpha_setbits = hashing::unique_in_range(0, nset, 0, nsite, true);
    auto beta_setbits = hashing::unique_in_range(1, nset, 0, nsite, true);
    mbf = {alpha_setbits, beta_setbits};
    ASSERT_EQ(mbf.nsetbit(), 2*nset);

    {
        /*
         * check that alpha channel iterator returns the correct site indices
         */
        auto it = alpha_setbits.cbegin();
        auto fn = [&it](size_t isite){
            ASSERT_EQ(isite, *(it++));
        };
        mbf.foreach_alpha(fn);
        ASSERT_EQ(it, alpha_setbits.cend());
    }
    {
        /*
         * check that beta channel iterator returns the correct site indices
         */
        auto it = beta_setbits.cbegin();
        auto fn = [&it](size_t isite){
            ASSERT_EQ(isite, *(it++));
        };
        mbf.foreach_beta(fn);
        ASSERT_EQ(it, beta_setbits.cend());
    }

#if 0
    defs::inds isites_openshell;
    std::cout << mbf << std::endl;
    for (size_t isite=0ul; isite<nsite; ++isite){
        if (mbf.get({0, isite})!=mbf.get({1, isite})) isites_openshell.push_back(isite);
    }
    //auto it = isites_openshell.cbegin();
    auto fn = [&](size_t isite){
        std::cout << isite << std::endl;
        //ASSERT_EQ(isite, (*it++));
    };
    std::cout << isites_openshell << std::endl;
    mbf.foreach_openshell(fn);
#endif
}