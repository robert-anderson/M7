//
// Created by Robert J. Anderson on 09/03/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/basis/Planewaves.h"

TEST(Planewaves, HeterogeneousExtents) {
    defs::uintv_t shape = {3, 5, 2};
    Planewaves planewaves(shape);
    size_t chk_size = 1;
    for (auto &extent: shape) chk_size*=(2*extent+1);
    ASSERT_EQ(planewaves.size(), chk_size);
    for (size_t i=0ul; i<planewaves.size(); ++i){
        auto& momvec = planewaves[i];
        ASSERT_EQ(planewaves.encode(momvec), i);
    }
}

TEST(Planewaves, HomogeneousExtents) {
    const size_t ndim=3, nwave = 4;
    Planewaves planewaves(ndim, nwave);
    ASSERT_EQ(planewaves.size(), std::pow(2*nwave+1, ndim));
    for (size_t i=0ul; i<planewaves.size(); ++i){
        auto& momvec = planewaves[i];
        ASSERT_EQ(planewaves.encode(momvec), i);
    }
}