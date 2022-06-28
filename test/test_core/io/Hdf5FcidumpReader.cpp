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
    auto shape = hdf5::DatasetReader::get_shape<hsize_t>(fr, "FOCK_INDEX");
    ASSERT_EQ(shape[0], 6ul);
    ASSERT_EQ(shape[1], 2ul);
    auto inds = fr.read_data<std::vector<int64_t>>("FOCK_INDEX");

    std::cout << inds << std::endl;

    hdf5::FileWriter fw("rja.h5");
    fw.write_data("FOCK_INDEX", inds);
}
