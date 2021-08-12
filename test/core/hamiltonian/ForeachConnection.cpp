//
// Created by rja on 06/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/hamiltonian/ForeachConnection.h"

TEST(ForeachConnection, FrmHubbard1d) {
    /*
     * loop over all connections using the general foreach class, and again with the Hubbard-optimized, and ensure that
     * the non-zero H-element connections returned are equivalent
     */
    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);
    conn::FrmOnv conn(ham.nsite());
    buffered::FrmOnv onv(ham.nsite());
    onv = {{0, 2, 4}, {0, 1, 5}};
    const defs::ham_t t = -1.0;

    std::vector<std::pair<defs::inds, defs::inds>> conns_general;
    foreach_conn::frm::Fermion foreach_general(ham);
    auto fn_general = [&](const conn::FrmOnv& conn, defs::ham_t helement){
        conns_general.emplace_back(conn.ann(), conn.cre());
        ASSERT_EQ(helement, t);
    };
    foreach_general.foreach<field::FrmOnv>(onv, fn_general, true);

    std::vector<std::pair<defs::inds, defs::inds>> conns_opt;
    foreach_conn::frm::Hubbard1D foreach_sym(ham);
    auto fn_opt = [&](const conn::FrmOnv& conn, defs::ham_t helement){
        conns_opt.emplace_back(conn.ann(), conn.cre());
        ASSERT_EQ(helement, t);
    };
    foreach_sym.foreach<field::FrmOnv>(onv, fn_opt, true);
    ASSERT_EQ(conns_opt, conns_general);
}


TEST(ForeachConnection, FrmNoSym) {
    Hamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);
    conn::FrmOnv conn(ham.nsite());
    buffered::FrmOnv onv(ham.nsite());
    ham.set_hf_mbf(onv, 0);

    auto counter = 0ul;
    auto fn = [&](const conn::FrmOnv& conn) {
        ++counter;
    };
    foreach_conn::frm::Fermion foreach(ham);

    const auto nocc = ham.nelec();
    const auto nvac = 2 * ham.nsite() - nocc;
    const auto nsingle = nocc * nvac;
    const auto ndouble = integer_utils::combinatorial(nocc, 2) * integer_utils::combinatorial(nvac, 2);

    foreach.foreach<field::FrmOnv>(onv, fn, false);
    ASSERT_EQ(counter, nsingle + ndouble);

    counter = 0ul;
    foreach.foreach<field::FrmOnv>(onv, fn, true);
    // can only excite to the vacant orbital to the right in both spin channels
    ASSERT_EQ(counter, 2ul);
}

#if 0
TEST(ForeachConnection, Hubbard1D) {
    FermionHamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);
    conn::Antisym<0> conn(ham.nsite());
    buffered::Onv<0> onv(ham.nsite());
    onv = {{0, 3, 5},
           {3, 4, 5}};
    /*
     * (100101,000111)
     * open boundary conditions:
     * alpha: 0->1, 3->2, 3->4, 5->4
     * beta:  3->2
     */
    auto counter = 0ul;
    auto fn_obc = [&](defs::ham_t helem) {
        ++counter;
    };
    foreach_conn::Hubbard1D foreach_obc(ham, conn, fn_obc, false);
    foreach_obc(onv);
    ASSERT_EQ(counter, 5ul);
    /*
     * with periodic boundary conditions:
     * beta 5->0
     */
    counter = 0ul;
    auto fn_pbc = [&](defs::ham_t helem) {
        ++counter;
    };
    foreach_conn::Hubbard1D foreach_pbc(ham, conn, fn_pbc, false);
    foreach_pbc(onv);
    ASSERT_EQ(counter, 5ul);
}


TEST(ForeachConnection, HubbardHolstein1D) {
    Hamiltonian<1> ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 4, 0.3, 1.0);
    conn::Antisym<1> conn(ham.nsite());
    buffered::Onv<1> onv(ham.nsite());
    onv.m_frm = {{0, 2},
                 {1, 3}};
    onv.m_bos = {0, 1, 4, 3};
    std::cout << onv << std::endl;
    /*
     * (1010,0101) [0 1 4 3 ]
     * open boundary conditions:
     * alpha: 0->1, 2->1, 2->3
     * beta:  1->0, 1->2, 3->2
     * bosons: (0 +1), (1 -1), (1 +1), (2 -1), (3 +1), (3 -1)
     * 3 + 3 + 6 = 12 total connections
     */
    auto counter = 0ul;
    auto fn_obc = [&](defs::ham_t helem) {
        ++counter;
    };
    foreach_conn::FermiBos foreach_obc(ham, conn,
                                       foreach_conn::Hubbard1D(ham, conn, fn_obc, false),
                                       fn_obc);
    foreach_obc(onv);
    ASSERT_EQ(counter, 12ul);
//    /*
//     * with periodic boundary conditions:
//     * beta 5->0
//     */
//    counter = 0ul;
//    auto fn_pbc = [&](defs::ham_t helem){
//        ++counter;
//    };
//    foreach_conn::Hubbard1D foreach_pbc(ham, conn, fn_pbc, false);
//    foreach_pbc(onv);
//    ASSERT_EQ(counter, 5ul);
}
#endif