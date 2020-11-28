//
// Created by rja on 19/07/2020.
//

#include <src/defs.h>
#include <src/core/io/FileReader.h>
#include "gtest/gtest.h"

TEST(FileReader, Skip){
    FileReader file_reader(defs::assets_root+"/RHF_N2_CCPVTZ/FCIDUMP");
    file_reader.skip(1);
}
