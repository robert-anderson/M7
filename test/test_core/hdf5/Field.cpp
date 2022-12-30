//
// Created by rja on 23/12/22.
//

#include "gtest/gtest.h"
#include "M7_lib/hdf5/IoManager.h"
#include "M7_lib/util/Hash.h"
#include "M7_lib/hdf5/File.h"
#include "M7_lib/wavefunction/FciInitializer.h"

TEST(FieldDataset, FrmOnvField) {
    HeisenbergFrmHam frm_ham(1.0, lattice::make("ortho", {10}, {1}));
    Hamiltonian ham(&frm_ham);
    // generate all spin determinants
    FciInitializer init(ham);
    {
        hdf5::FileWriter fw("tmp.h5");
        init.m_mbf_order_table.save(fw, "mbf_table", mpi::i_am_root());
    }
    FciInitializer::mbf_order_table_t table("mbf_table_load", init.m_mbf_order_table.m_row);
//    {
//        hdf5::FileReader fr("tmp.h5");
//        table.load(fr, );
//    }
}