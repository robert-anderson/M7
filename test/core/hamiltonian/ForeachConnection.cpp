//
// Created by rja on 06/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/hamiltonian/ForeachConnection.h"

namespace foreach_connection_test {
    typedef SingleFieldRow<fields::Onv<1>> result_row_t;
    typedef BufferedTable<result_row_t, 1> result_table_t;

    result_table_t hubbard_6site() {
        result_table_t table("", {{6}, 10});
        auto& row = table.m_row;
        row.push_back_jump();
        row.m_field = {{1, 4, 5}, {1, 3, 0, 0, 2, 1}};
        return table;
    }
}

TEST(ForeachConnection, FermionNoSym){
    FermionHamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);
    conn::Antisym<0> conn(ham.nsite());
    buffered::Onv<0> onv(ham.nsite());
    ham.set_hf_onv(onv, 0);

    auto counter = 0ul;
    auto fn = [&](defs::ham_t helem){
        ++counter;
    };
    foreach_conn::Fermion nosym(ham, conn, fn, true, false);
    nosym(onv);
    const auto nocc = ham.nelec();
    const auto nvac = 2*ham.nsite()-nocc;
    const auto nsingle = nocc*nvac;
    const auto ndouble = integer_utils::combinatorial(nocc, 2)*integer_utils::combinatorial(nvac, 2);
    ASSERT_EQ(counter, nsingle+ndouble);
    counter = 0ul;
    auto fn_hnonzero = [&](defs::ham_t helem){
        ++counter;
        ASSERT_FLOAT_EQ(helem, -1.0);
    };
    foreach_conn::Fermion nosym_hnonzero(ham, conn, fn_hnonzero);
    nosym_hnonzero(onv);
    ASSERT_EQ(counter, 2ul);
}


TEST(ForeachConnection, FermionHubbard1D){


#if 0
    FermionHamiltonian ham(defs::assets_root + "/Hubbard_U4_6site/FCIDUMP", 0);
    conn::Antisym<0> conn(ham.nsite());
    buffered::Onv<0> onv(ham.nsite());
    ham.set_hf_onv(onv, 0);

    auto counter = 0ul;
    auto fn = [&](defs::ham_t helem){
        ++counter;
    };
    foreach_conn::Fermion nosym(ham, conn, fn, true, false);
    nosym(onv);
    const auto nocc = ham.nelec();
    const auto nvac = 2*ham.nsite()-nocc;
    const auto nsingle = nocc*nvac;
    const auto ndouble = integer_utils::combinatorial(nocc, 2)*integer_utils::combinatorial(nvac, 2);
    ASSERT_EQ(counter, nsingle+ndouble);
    counter = 0ul;
    auto fn_hnonzero = [&](defs::ham_t helem){
        ++counter;
        ASSERT_FLOAT_EQ(helem, -1.0);
    };
    foreach_conn::Fermion nosym_hnonzero(ham, conn, fn_hnonzero);
    nosym_hnonzero(onv);
    ASSERT_EQ(counter, 2ul);
#endif
}