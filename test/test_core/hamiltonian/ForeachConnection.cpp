//
// Created by Robert J. Anderson on 06/06/2021.
//

#include <M7_lib/foreach/Foreach.h>
#include <M7_lib/caches/DecodedMbf.h>
#include "gtest/gtest.h"


//struct ForEachInExsig {
//    const size_t m_exsig;
//
//    virtual void loop() = 0;
//
//    virtual bool draw(const field::FrmOnv &src_onv,
//                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
//                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);
//};



//struct ForeachGroup {
//    typedef std::function<void(conn::FrmOnv)> body_fn_t;
//    void loop()
//};


#if 0
TEST(ForeachConnection, FrmNoSymDoubles){
    Hamiltonian ham("/home/rja/CLionProjects/M7/assets/HF_RDMs/FCIDUMP", 0);
    FrmConserveExcitIter iter(exsig_utils::ex_double, ham);

    auto fn = [](const field::FrmOnv& dst, defs::ham_t h){
        std::cout << dst << " " << h << std::endl;
    };
    buffered::FrmOnv onv(ham.nsite());
    onv = {{0, 1, 2}, {0, 1, 2}};

    iter.foreach<field::FrmOnv>(onv, fn, true);
}

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


TEST(ForeachConnection, HubbardHolstein) {
    auto fname = defs::assets_root + "/Hubbard_U4_3site/FCIDUMP";
    auto fname_eb = defs::assets_root + "/Hubbard_U4_3site/EBDUMP_HH_V1.4";
    auto fname_bos = defs::assets_root + "/Hubbard_U4_3site/BOSDUMP_HH_W0.3";
    Hamiltonian ham(fname, fname_eb, fname_bos, false, 3);

    conn::FrmBosOnv conn(ham.nsite());
    buffered::FrmBosOnv onv(ham.nsite());

    onv = {{1, 2, 3, 4}, {1, 0, 2}};

    struct Result {
        defs::inds m_frm_cre, m_frm_ann, m_bos_cre, m_bos_ann;
    };

    std::vector<Result> results;
    auto fn = [&](const conn::FrmBosOnv& conn) {
        results.push_back({conn.m_frm.m_cre.inds(), conn.m_frm.m_ann.inds(),
                             conn.m_bos.m_cre.to_vector(), conn.m_bos.m_ann.to_vector()});
        std::cout << results.back().m_bos_ann << " " << results.back().m_bos_cre << std::endl;
    };
    foreach_conn::frm_bos::Holstein foreach(ham);

    foreach.foreach<field::FrmBosOnv>(onv, fn, false);
}
#endif


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