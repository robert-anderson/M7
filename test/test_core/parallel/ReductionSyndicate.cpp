//
// Created by Robert J. Anderson on 19/11/2020.
//

#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/parallel/Reduction.h>
#include "gtest/gtest.h"
#include "M7_lib/parallel/ReductionMember.h"

/*
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
 */


TEST(ReductionSyndicate, New){
    //NdReduction<int, 2> red({2, 3});

}