//
// Created by rja on 26/01/23.
//

#include "test_core/defs.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/table/GlobalAccumulation.h"

TEST(GlobalAccumulation, Test) {
    const uint_t nchar = 10ul;
    const uint_t nrank_pair = (mpi::nrank()*(mpi::nrank()-1))/2;
    struct row_t : Row {
        field::String m_key;
        field::Number<int> m_value;
        row_t(): m_key(this, nchar), m_value(this){}
        field::String &key_field() {
            return m_key;
        };
        field::Number<int> &value_field() {
            return m_value;
        };
    };
    GlobalAccumulation<row_t> ga("global accumulator test", row_t{});
    ASSERT_EQ(ga.naccum(), 0ul);
    buffered::String work_key(nchar);
    buffered::Number<int> work_value;
    work_key = "abc";
    work_value = 4 * mpi::irank() + 23;
    ga.add(work_key, work_value);
    work_key = "def";
    work_value = -3 * mpi::irank() + -56;
    ga.add(work_key, work_value);
    ga.update();
    ASSERT_EQ(ga.naccum(), 1ul);
    auto& row = ga.current().m_row;
    row.jump(0);
    ASSERT_EQ(row.m_key, "abc");
    ASSERT_EQ(row.m_value[0], 4 * nrank_pair + 23 * mpi::nrank());
    row.jump(1);
    ASSERT_EQ(row.m_key, "def");
    ASSERT_EQ(row.m_value[0], -3 * nrank_pair + -56 * mpi::nrank());
}