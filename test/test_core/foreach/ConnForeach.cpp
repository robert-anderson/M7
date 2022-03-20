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

    auto fn = [](const conn::FrmOnv& conn){
        std::cout << conn.m_cre.inds() << " " << conn.m_ann.inds() << std::endl;
    };
    conn_foreach::frm::General<2> foreach(nsite, fn);
    foreach.loop(mbf);

    std::cout << mbf << std::endl;
}
//TEST(ConnForeach, FrmBosLatticeSingles){
//
//
//}