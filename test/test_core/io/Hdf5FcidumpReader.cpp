//
// Created by anderson on 24/06/2022.
//


#include "test_core/defs.h"
#include "M7_lib/io/FcidumpTextFileReader.h"
#include "M7_lib/hdf5/DatasetTransaction.h"

TEST(Hdf5FcidumpReader, Header) {
    using namespace hdf5;
    hdf5::FileReader fr(PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5");
    auto shape = DatasetLoader::read_format(fr, "FOCK_INDEX", false, true).m_h5_shape;
    ASSERT_EQ(shape[0], 6ul);
    ASSERT_EQ(shape[1], 2ul);
    auto inds = DatasetLoader::load_vector<int64_t>(fr, "FOCK_INDEX", DatasetLoader::Options());

    hdf5::FileWriter fw("rja.h5");
    hdf5::DatasetSaver::save_vector(fw, "FOCK_INDEX", inds, 0ul);
}
