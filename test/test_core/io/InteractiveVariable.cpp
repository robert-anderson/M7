//
// Created by Robert J. Anderson on 16/03/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/io/InteractiveVariable.h"

TEST(InteractiveVariable, Test){
    InteractiveVariableFile iv("some_uint");
    std::vector<float> i(4);
    iv.read(i);
}