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

    hdf5::FileWriter fw("tmp.h5");
    init.m_mbf_order_table.m_row.m_field.save(fw, "frm_onv", mpi::i_am_root());
}