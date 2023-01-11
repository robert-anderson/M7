//
// Created by rja on 23/12/22.
//

#include "gtest/gtest.h"
#include "M7_lib/hdf5/DatasetFormat.h"
#include "M7_lib/util/Hash.h"
#include "M7_lib/hdf5/File.h"
#include "M7_lib/wavefunction/FciInitializer.h"

TEST(FieldDataset, SaveEmpty) {
    typedef SingleFieldRow<field::Number<int>> row_t;
    buffered::Table<row_t> table("numbers");
    ASSERT_EQ(table.nrow_in_use(), 0ul);
    {
        hdf5::FileWriter fw("tmp.h5");
        table.save(fw, "empty_table", mpi::i_am_root());
    }
}

TEST(FieldDataset, FrmOnvField) {
    HeisenbergFrmHam frm_ham(1.0, lattice::make("ortho", {10}, {1}));
    Hamiltonian ham(&frm_ham);
    // generate all spin determinants
    FciInitializer init(ham);
    ASSERT_EQ(init.m_mbf_order_table.nrow_in_use(), 252ul);
    const uint_t max_nitem_per_op = 12ul;
    {
        hdf5::FileWriter fw("tmp.h5");
        init.m_mbf_order_table.save(fw, "mbf_table", max_nitem_per_op, mpi::i_am_root());
    }
    FciInitializer::mbf_order_table_t table("mbf_table_load", init.m_mbf_order_table.m_row);
    {
        hdf5::FileReader fr("tmp.h5");
        table.load(fr, "mbf_table", false, true);
    }
    // check that the MappedTable has been correctly restored
    auto row = table.m_row;
    for (row.restart(); row; ++row){
        auto res = table.lookup(row.m_field);
        ASSERT_TRUE(res);
        ASSERT_EQ(row.index(), res.index());
    }
}