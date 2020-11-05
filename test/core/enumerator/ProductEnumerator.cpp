//
// Created by rja on 07/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/enumerator/ProductEnumerator.h"
#include "src/core/util/utils.h"

TEST(ProductEnumerator, Test){
    const size_t nind = 4;
    const size_t extent = 5;

    ProductEnumerator enumerator(nind, extent);
    defs::inds inds(nind);
    size_t i=~0ul;
    while(enumerator.next(inds, i)){}
    ASSERT_EQ(i, std::pow(extent, nind));
}