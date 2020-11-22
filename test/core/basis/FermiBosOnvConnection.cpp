//
// Created by RJA on 21/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/field/Elements.h"

TEST(FermiBosConnection, Test){
    const elements::FermiBosOnv fonv(4);
    fonv.print();
}