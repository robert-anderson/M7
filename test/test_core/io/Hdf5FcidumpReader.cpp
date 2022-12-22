//
// Created by anderson on 24/06/2022.
//

#include <utility>
#include <M7_lib/hdf5/Dataset.h>

#include "test_core/defs.h"
#include "M7_lib/io/FcidumpTextFileReader.h"

using namespace hdf5;

TEST(Hdf5FcidumpReader, Header) {
    hdf5::FileReader fr(PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5");
    auto shape = hdf5::DatasetReader::get_hdf5_shape(fr, "FOCK_INDEX");
    ASSERT_EQ(shape[0], 6ul);
    ASSERT_EQ(shape[1], 2ul);
    auto inds = fr.load_dataset<v_t<int64_t>>("FOCK_INDEX", false, true);

    hdf5::FileWriter fw("rja.h5");
    fw.save_dataset("FOCK_INDEX", inds, mpi::i_am_root());
}
