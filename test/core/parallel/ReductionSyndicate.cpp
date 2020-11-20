//
// Created by RJA on 19/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/parallel/ReductionMember.h"

TEST(ReductionSyndicate, Test){
    ReductionSyndicate syndicate;
    ReductionMember<double, 1> member(syndicate, {2});
    member(1) = 2;
    std::cout << member(1) << std::endl;
    std::cout << member.reduced(1) << std::endl;
    syndicate.sum();
    std::cout << member.reduced(1) << std::endl;
    syndicate.zero();
    std::cout << member.reduced(1) << std::endl;
}