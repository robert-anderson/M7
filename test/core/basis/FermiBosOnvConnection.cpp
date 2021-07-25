//
// Created by RJA on 21/11/2020.
//

#include <src/core/table/BufferedFields.h>
#include "gtest/gtest.h"

TEST(FermiBosConnection, Test){
    const buffered::FrmBosOnv fonv(4);
    std::cout << fonv.to_string() << std::endl;
}