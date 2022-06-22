//
// Created by Robert J. Anderson on 12/4/21.
//

#include "gtest/gtest.h"
#include "M7_lib/io/FortranNamelistReader.h"

TEST(FortranNamelistReader, IsolateValue) {
    std::string line = " &FCI NORB=  12,NELEC=12,MS2=0, TREL=.TRUE., ORBSYM=3,2,1,1,4,1,5,3,2,1,5,6";
    ASSERT_EQ(FortranNamelistReader::isolate_value(line, "NORB"), "12");
    ASSERT_EQ(FortranNamelistReader::isolate_value(line, "NELEC"), "12");
    ASSERT_EQ(FortranNamelistReader::isolate_value(line, "TREL"), ".TRUE.");
    ASSERT_EQ(FortranNamelistReader::isolate_value(line, "ORBSYM"), "3,2,1,1,4,1,5,3,2,1,5,6");
}

TEST(FortranNamelistReader, FromFile) {
    FortranNamelistReader header(PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP");
    /*
     * &FCI NORB=   6,NELEC= 6,MS2=0,
     * ORBSYM=1,3,2,6,7,5
     * ISYM=1,
     * &END
     */
    ASSERT_EQ(header.read_uint("NORB"), 6);
    ASSERT_EQ(header.read_uint("NELEC"), 6);
    ASSERT_EQ(header.read_int("MS2", 2), 0);
    defs::inds_t orbsym = {0, 2, 1, 5, 6, 4};
    ASSERT_EQ(header.read_uints("ORBSYM", -1), orbsym);
}