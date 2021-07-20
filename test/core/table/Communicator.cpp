//
// Created by rja on 18/05/2021.
//

#include <src/core/table/BufferedFields.h>
#include "src/core/table/Communicator.h"
#include "gtest/gtest.h"

namespace communicator_test {
    struct TestRow : Row {
        fields::Number<size_t> m_key;
        fields::Number<double> m_value;
        TestRow(): m_key(this, "key"), m_value(this, "value"){}
        fields::Number<size_t> &key_field() {
            return m_key;
        };
    };
}
TEST(Communicator, SharedRow) {
    using namespace communicator_test;
    typedef Communicator<TestRow, TestRow> comm_t;
    comm_t comm("test communicator", 0.5, 5, 1, {{}, 10}, {{}}, 0.04, 20ul);
    const size_t nrow_per_rank_expect = 6;
    auto &row = comm.m_store.m_row;
    for (size_t i = 0; i<nrow_per_rank_expect*mpi::nrank(); ++i){
        auto key = 123+i*5;
        if (!mpi::i_am(comm.m_ra.get_rank(key))) continue;
        row.push_back_jump();
        row.m_key = key;
        row.m_value = 2.8*i;
    }
    ASSERT_EQ(mpi::all_sum(comm.m_store.m_hwm), nrow_per_rank_expect*mpi::nrank());

    comm_t::SharedRow shared_row(comm, {0, 0}, "test shared row");
}