//
// Created by anderson on 03/06/2022.
//

#ifndef M7_TEST_DEFS_H
#define M7_TEST_DEFS_H

#include "gtest/gtest.h"

#define ASSERT_COMPLEX_EQ(v1, v2) ASSERT_FLOAT_EQ(std::abs(v1-v2), 0.0);
#define ASSERT_NEARLY_EQ(v1, v2, tol) ASSERT_LE(std::abs(v1-v2), std::abs(tol));

#endif //M7_TEST_DEFS_H