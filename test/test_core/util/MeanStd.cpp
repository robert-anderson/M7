//
// Created by rja on 23/06/22.
//

#include "test_core/defs.h"
#include "M7_lib/util/MeanStd.h"

TEST(UtilMeanStd, MeanAndStd) {
    MeanStd<double> mean_std({1, 2, 3.8, 4, -0.35, 0.6});
    ASSERT_NEAR_EQ(mean_std.m_mean, 1.8416666666666668);
    ASSERT_NEAR_EQ(mean_std.m_std, 1.6110081384717527);
}
