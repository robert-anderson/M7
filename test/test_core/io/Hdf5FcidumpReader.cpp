//
// Created by anderson on 24/06/2022.
//


#include "test_core/defs.h"
#include "M7_lib/io/FcidumpTextFileReader.h"
#include "M7_lib/hdf5/DatasetSaver.h"

using namespace hdf5;

TEST(Hdf5FcidumpReader, Header) {
    hdf5::FileReader fr(PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5");
    auto shape = fr.get_dataset_shape("FOCK_INDEX");
    ASSERT_EQ(shape[0], 6ul);
    ASSERT_EQ(shape[1], 2ul);
    auto inds = fr.load_dataset<v_t<int64_t>>("FOCK_INDEX", false, true);

    hdf5::FileWriter fw("rja.h5");
    hdf5::DatasetSaver::save_vector(fw, "FOCK_INDEX", inds, 0ul);
}
