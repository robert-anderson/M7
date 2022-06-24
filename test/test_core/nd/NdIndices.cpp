//
// Created by Robert J. Anderson on 19/08/2021.
//

#include "M7_lib/nd/NdIndices.h"
#include "gtest/gtest.h"

TEST(NdIndices, Scalar){
    NdIndices<0> nd_inds;
    ASSERT_EQ(nd_inds.size(), 1ul);
}

TEST(NdIndices, Vector){
    const uint_t n = 12;
    NdIndices<1> nd_inds(n);
    for (uint_t i=0ul; i<n; ++i)
        ASSERT_EQ(nd_inds[i][0], i);
}

TEST(NdIndices, MultiDim){
    const uinta_t<4> shape = {4, 3, 5, 7};
    uinta_t<4> inds{};
    NdIndices<4> nd_inds(shape);
    uint_t iflat = 0ul;
    for (inds[0]=0ul; inds[0]<shape[0]; ++inds[0])
        for (inds[1]=0ul; inds[1]<shape[1]; ++inds[1])
            for (inds[2]=0ul; inds[2]<shape[2]; ++inds[2])
                for (inds[3]=0ul; inds[3]<shape[3]; ++inds[3])
                    ASSERT_EQ(nd_inds[iflat++], inds);
}