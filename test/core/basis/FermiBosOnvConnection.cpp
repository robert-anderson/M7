//
// Created by RJA on 21/11/2020.
//

#include <src/core/field/BufferedFields.h>
#include "gtest/gtest.h"

TEST(FermiBosConnection, Test){
    const buffered::FermiBosOnv fonv(4);
    std::cout << fonv.to_string() << std::endl;
}