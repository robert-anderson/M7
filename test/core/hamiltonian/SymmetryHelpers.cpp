//
// Created by rja on 01/06/2021.
//

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/hamiltonian/ForeachConnection.h>
#include "gtest/gtest.h"
#include "src/core/hamiltonian/SymmetryHelpers.h"

/*
 * foreach loops over connections can be very useful and reusable.
 *
 * the most flexible implementation can be prototyped as follows:
 *
 * // the hamiltonian object
 * const Hamiltonian<enable_bosons>& ham
 * // a reference to the source ONV from which connections are to be generated. this could be for either ONV type
 * const fields::Onv<enable_bosons>& src_onv
 * // a connection object corresponding to the ONV type of src_onv. this connection is the receptacle for all generated
 * // connections from src_onv, so we need to be able to access it by reference from both the looping scope, and the
 * // closure
 * conn::Antisym<enable_bosons> conn;
 * // a capturing lambda defining what to do with each connection generated:
 * auto body_fn = [&](const defs::ham_t& helement){
 *      // do whatever with conn, which has been set to a valid connection of the hamiltonian
 * }
 *
 * // the foreach function sets conn state to that of next connection and then calls body_fn
 * foreach(ham, src_onv, conn, body_fn);
 *
 */


TEST(SymmetryHelpers, Test) {
    Hamiltonian<1> ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0, 1, 1.4, 0.3);
    ham_sym_helpers::FermiBos sym_helper(ham);
    BufferedTable<SingleFieldRow<fields::Onv<1>>> table("", {{ham.nsite()}});
    buffered::Onv<1> src_onv(ham.nsite());
    src_onv.m_frm = {0, 3, 5, 7};
    src_onv.m_bos = {1, 0, 2, 3};
    auto body_fn = [&](const conn::Antisym<1> & conn, const fields::Onv<1> &dst_onv, const defs::ham_t &helem){
        table.m_row.push_back_jump();
        table.m_row.m_field = dst_onv;
    };
    sym_helper.foreach_connection(src_onv, body_fn, false, false, false);
}


TEST(SymmetryHelpers, Test2) {

    Hamiltonian<0> ham(defs::assets_root + "/Hubbard_U4_4site/FCIDUMP", 0);
    buffered::Onv<0> src_onv(ham.nsite());
    ham.set_hf_onv(src_onv, 0);
    conn::Antisym<0> conn(ham.nsite());

    buffered::Onv<0> dst_onv(ham.nsite());
    auto body_fn = [&](defs::ham_t helement) {
        conn.apply(src_onv, dst_onv);
    };
    foreach_conn::Fermion(ham, conn, body_fn, true, false)(src_onv);

}