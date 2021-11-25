//
// Created by rja on 25/11/2021.
//

#include "gtest/gtest.h"
#include "src/core/wavefunction/WalkerTable.h"

TEST(WalkerTable, Fields){
    WalkerTable table(WalkerTableRow({0, 5}, 1, 1, false));
    auto& row = table.m_row;
    std::vector<const FieldBase*> fields = {
            &row.m_mbf, &row.m_weight, &row.m_hdiag, &row.m_initiator, &row.m_deterministic, &row.m_ref_conn};
    size_t i = 0ul;
    for (auto field: fields) {
        ASSERT_EQ(field->m_row, &row);
        ASSERT_EQ(field, row.m_fields[i++]);
    }
}