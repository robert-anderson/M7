//
// Created by Robert J. Anderson on 19/07/2020.
//

#include <M7_lib/defs.h>
#include <M7_lib/io/FileReader.h>
#include "gtest/gtest.h"

TEST(FileReader, Skip){
    FileReader file_reader(PROJECT_ROOT"/assets/RHF_N2_CCPVTZ/FCIDUMP");
    uint_t n = 0ul;
    while (file_reader.next()) ++n;
    ASSERT_EQ(n, 180704);
}
