//
// Created by anderson on 17/02/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/foreach/ConnForeach.h"

TEST(ConnForeach, FrmGeneralEx1100){
    const size_t nsite = 8;
    defs::inds setbits = {1, 4, 6, 9, 12};
    buffered::FrmOnv mbf(nsite);
    mbf = setbits;
    defs::inds clrbits;
    for (size_t ibit = 0ul; ibit<2*nsite; ++ibit) if (!mbf.get(ibit)) clrbits.push_back(ibit);

    size_t iiter = 0ul;
    auto fn = [&](const conn::FrmOnv& conn){
        const auto cre = conn.m_cre[0];
        const auto ann = conn.m_ann[0];
        const auto iann = iiter/clrbits.size();
        const auto icre = iiter-iann*clrbits.size();
        ASSERT_EQ(cre, clrbits[icre]);
        ASSERT_EQ(ann, setbits[iann]);
        ++iiter;
    };
    conn_foreach::frm::General<1> foreach(nsite);
    ASSERT_EQ(foreach.m_exsig, exsig_utils::ex_single);
    foreach.loop_fn(mbf, fn);
    ASSERT_EQ(iiter, setbits.size()*clrbits.size());
}