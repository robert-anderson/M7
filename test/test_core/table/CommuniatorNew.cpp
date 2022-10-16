//
// Created by Robert J. Anderson on 18/05/2021.
//

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/communication/Communicator.h"
#include "gtest/gtest.h"

namespace communicator_new_test {
    struct TestRow : Row {
        field::Number<uint_t> m_key;
        field::Number<double> m_value;
        TestRow(): m_key(this, "key"), m_value(this, "value"){}
        field::Number<uint_t> &key_field() {
            return m_key;
        };
    };
}

TEST(CommunicatorNew, SharedRow) {
    using namespace communicator_new_test;
    typedef communicator::BasicSend<TestRow, TestRow> comm_t;
    const Sizing store_sizing = {20, 0.0};
    const Sizing comm_sizing = {20, 0.0};
    const DistribOptions dist_opts;
    comm_t comm("test communicator", {}, dist_opts, store_sizing, {}, comm_sizing);
    const uint_t nrow_per_rank_expect = 6;
    auto &row = comm.m_store.m_row;
    for (uint_t i = 0; i<nrow_per_rank_expect*mpi::nrank(); ++i){
        BufferedField<field::Number<uint_t>> key;
        key = 123+i*5;
        if (!mpi::i_am(comm.m_dist.irank(key))) continue;
        row.push_back_jump();
        row.m_key = key;
        row.m_value = 2.8*i;
    }
    ASSERT_EQ(mpi::all_sum(comm.m_store.nrow_in_use()), nrow_per_rank_expect*mpi::nrank());

//    comm_t::SharedRow shared_row(comm, {0, 0}, "test shared row");
}

TEST(CommunicatorNew, Redistribution) {
    using namespace communicator_new_test;
    typedef communicator::BasicSend<TestRow, TestRow> comm_t;
    const Sizing store_sizing = {20, 0.0};
    const Sizing comm_sizing = {20, 0.0};
    const DistribOptions dist_opts;
    comm_t comm("test communicator", {}, dist_opts, store_sizing, {}, comm_sizing);
    const uint_t nrow_per_rank_expect = 7;
    auto &row = comm.m_store.m_row;
    v_t<double> work_figures(comm.m_dist.nblock());
    for (uint_t i = 0; i<nrow_per_rank_expect*mpi::nrank(); ++i){
        BufferedField<field::Number<uint_t>> key;
        key = 123+i*5;
        if (!mpi::i_am(comm.m_dist.irank(key))) continue;
        comm.m_store.insert(key);
        row.m_key = key;
        row.m_value = 2.8*i;
        if (!(i%4)) row.protect();
        if (!(i%2)) row.protect();
        /*
         * make up an amount of work done
         */
        double work_done = hash::in_range(i, 2, 100)/double(100);
        comm.m_store.accumulate_work_figure(key, work_done);
        work_figures[comm.m_dist.iblock(key)]+=work_done;
    }
    //Distribution orig_dist(work_figures.size(), mpi::nrank());
    //Redistributor redist_chk(orig_dist.block_iranks(), work_figures, mpi::nrank());

    comm.m_store.redistribute();

    for (uint_t i = 0; i<nrow_per_rank_expect*mpi::nrank(); ++i){
        BufferedField<field::Number<uint_t>> key;
        key = 123+i*5;
        if (!mpi::i_am(comm.m_dist.irank(key))) continue;
        comm.m_store.lookup(key, row);
        if (!(i%4)) {
            ASSERT_EQ(row.protection_level(), 2ul);
        }
        else if (!(i%2)) {
            ASSERT_EQ(row.protection_level(), 1ul);
        }
    }

//    comm_t::SharedRow shared_row(comm, {0, 0}, "test shared row");
}