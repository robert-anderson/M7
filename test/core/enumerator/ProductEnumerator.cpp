//
// Created by rja on 07/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/enumerator/ProductEnumerator.h"
#include "src/core/util/utils.h"

TEST(ProductEnumerator, Test){
    const size_t extent = 5;
    ProductEnumerator<4> enumerator(extent);
    defs::inds inds(4);
    size_t iflat=0ul;
    while(enumerator.next(inds)){
        iflat++;
    }
    ASSERT_EQ(iflat, std::pow(extent, 4));
}