//
// Created by rja on 19/07/2020.
//

#include "gtest/gtest.h"
#include "src/core/parallel/Reducible.h"

TEST(Reducible, MaxLocMinLoc){
    typedef double T;
    Reducible<T> reducible;
    auto to_T = [](size_t irank){
        return std::pow(1.13, (irank+7*3)%5+1);
    };
    reducible = to_T(mpi::irank());
    std::pair<T, int> max{std::numeric_limits<T>::min(), 0};
    std::pair<T, int> min{std::numeric_limits<T>::max(), 0};
    for (size_t i=0ul; i<mpi::nrank(); ++i) {
        auto d = to_T(i);
        if (d>max.first) max = {d, i};
        if (d<min.first) min = {d, i};
    }

    reducible.mpi_maxloc();
    ASSERT_EQ(max.first, reducible.reduced());
    ASSERT_EQ(max.second, reducible.irank());

    reducible.mpi_minloc();
    ASSERT_EQ(min.first, reducible.reduced());
    ASSERT_EQ(min.second, reducible.irank());

}