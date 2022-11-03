//
// Created by rja on 02/11/22.
//

#include "test_core/defs.h"

TEST(FpTol, NearEqual){
    const double b = 0.5;
    const double rtol = 0.03125;
    const double atol = 0.0078125;
    ASSERT_TRUE(fptol::near_equal(b, b, rtol, atol));
    // near_equal condition:
    // std::abs(a-b) <= (atol+rtol*std::abs(b));
    // i.e. a <= atol + b + rtol*|b|    if a > b
    //      a >= b - atol - rtol*|b|    else
    const auto ahi = atol + b + rtol*std::abs(b);
    const auto alo = b - atol - rtol*std::abs(b);
    const auto delta = 0.0000000001;//std::numeric_limits<double>::epsilon();

    ASSERT_TRUE(fptol::near_equal(alo, b, rtol, atol));
    ASSERT_TRUE(fptol::near_equal(ahi, b, rtol, atol));

    ASSERT_TRUE(fptol::near_equal(alo + delta, b, rtol, atol));
    ASSERT_TRUE(fptol::near_equal(ahi - delta, b, rtol, atol));

    ASSERT_FALSE(fptol::near_equal(alo - delta, b, rtol, atol));
    ASSERT_FALSE(fptol::near_equal(ahi + delta, b, rtol, atol));
}

TEST(FpTol, NearZero){
    const double b = std::pow(2.0, -11.0);
    const double tol = std::pow(2.0, -10.0);
    ASSERT_TRUE(fptol::near_zero(b, tol));
    ASSERT_TRUE(fptol::near_zero(-b, tol));
    ASSERT_FALSE(fptol::near_zero(b+2*tol, tol));
    ASSERT_FALSE(fptol::near_zero(b-2*tol, tol));
}