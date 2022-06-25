//
// Created by anderson on 24/06/2022.
//

#include "test_core/defs.h"
#include "M7_lib/io/HDF5Wrapper.h"
#include "M7_lib/io/FcidumpFileReader.h"

TEST(Hdf5FcidumpReader, Header) {
    hdf5::FileReader fr(PROJECT_ROOT"/assets/N2_Molcas/molcas.FciDmp.h5");
    auto s = fr.read_attr<std::string>("MOLCAS_MODULE", std::string("123"));
    std::cout << s << "|" << std::endl;
//    FcidumpInfo info(fr);
}