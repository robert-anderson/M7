//
// Created by rja on 16/03/2021.
//

#include "gtest/gtest.h"
#include "src/core/io/InteractiveVariable.h"

TEST(InteractiveVariable, Test){
    InteractiveVariableFile iv("some_uint");
    std::vector<float> i(4);
    iv.read(i);
    utils::print(i);
}