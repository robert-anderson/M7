//
// Created by Robert J. Anderson on 25/11/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/wavefunction/WalkerTable.h"

TEST(WalkerTable, Fields){
    sys::Basis basis(5, 7);
    WalkerTable table(WalkerTableRow(basis, 1, 1, false));
    auto& row = table.m_row;

    /**
     * use overloading to handle the MultiField case of FrmBosOnv
     */
    struct Appender {
        static void append(std::vector<const FieldBase*>& v, const FieldBase& f){
            v.push_back(&f);
        }
        static void append(std::vector<const FieldBase*>& v, const field::FrmBosOnv& f){
            v.push_back(&f.m_frm);
            v.push_back(&f.m_bos);
        }
    };

    std::vector<const FieldBase*> fields;
    Appender::append(fields, row.m_mbf);
    fields.insert(fields.end(), {&row.m_weight, &row.m_hdiag, &row.m_initiator, &row.m_deterministic, &row.m_ref_conn});
    size_t i = 0ul;
    for (auto field: fields) {
        ASSERT_EQ(field->m_row, &row);
        ASSERT_EQ(field, row.m_fields[i++]);
    }
}
