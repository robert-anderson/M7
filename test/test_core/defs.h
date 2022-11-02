//
// Created by anderson on 03/06/2022.
//

#ifndef M7_TEST_DEFS_H
#define M7_TEST_DEFS_H

#include "gtest/gtest.h"
#include "M7_lib/util/FpTol.h"
#include "M7_lib/io/Logging.h"

namespace test_defs {
    template<typename T>
    struct NearlyEqHelper {
        const T m_a;
        const arith::comp_t<T> m_rtol;
        const arith::comp_t<T> m_atol;
        bool operator==(const T& b) const {
            return fptol::nearly_equal(m_a, b, m_rtol, m_atol);
        }

        /**
         * overload this operator so that the value is properly represented as a string by GTest upon ASSERT failure
         */
        friend std::ostream& operator<<(std::ostream& os, const NearlyEqHelper& v) {
            return os << v.m_a;
        }
    };
    /**
     * can't instantiate directly from ctor without specifying the type, so use a maker to work around this
     */
    template<typename T>
    NearlyEqHelper<T> make_helper(T a, arith::comp_t<T> rtol, arith::comp_t<T> atol) {
        return {a, rtol, atol};
    }
}

/*
 * GTest doesn't define a very complete means of comparing floats, and doesn't deal with complex types at all, so here
 * the ASSERT_EQ macro is used in conjunction with an instance of a helper class in order to perform the necessary
 * near or numeric comparison, but still allow Gtest to output the values upon failure
 */
#define ASSERT_NEARLY_EQ_TOL(a, b, rtol, atol) ASSERT_EQ(test_defs::make_helper(a, rtol, atol), b);
#define ASSERT_NEARLY_EQ(a, b) ASSERT_NEARLY_EQ_TOL(a, b, fptol::default_rtol(a), fptol::default_atol(a));
#define ASSERT_NEARLY_ZERO_TOL(b, atol) ASSERT_NEARLY_EQ_TOL(0.0, b, 0.0, atol);
#define ASSERT_NEARLY_ZERO(b) ASSERT_NEARLY_ZERO_TOL(b, fptol::default_ztol(b));

// TODO: make these work for vectors

#endif //M7_TEST_DEFS_H