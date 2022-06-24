//
// Created by Robert J. Anderson on 29/01/2021.
//

#ifndef M7_MPIASSERT_H
#define M7_MPIASSERT_H

#include "MPIWrapper.h"
#include <M7_lib/util/FpTol.h>
#include <M7_lib/io/Logging.h>

/**
 * This file furnishes a set of macros for code verification.
 * Most of these are defined for binary operators, in which cases the left- and right-hand side expressions must be
 * provided along with a reason why (lhs) OP (rhs) would evaluate to false.
 *
 * ASSERT macros only evaluate the expression (lhs) OP (rhs) when NDEBUG is undefined
 * REQUIRE macros are compiled in every build
 *
 * There is an additional unary pair of macros to test for the truth of a boolean, in this case the values of
 * and individual operands cannot be reported in the error log
 *
 * Macros with the COLLECTIVE_ prefix are called on all ranks, and if the expression evaluates to false on any rank,
 * MPI finalization will be called on all ranks. Macros with the INDEPENDENT_ prefix do not require any synchronization
 * and any rank which evaluates the expression to false will abort. COLLECTIVE_ macros are favoured where applicable.
 */
#define MPI_ABORT(msg) {log::error("MPI_ABORT_ called at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}


namespace asserts {

    static void abort(const char *file, int line, bool collective, const std::string &reason) {
        if (collective || mpi::nrank() == 1) {
            log::error("Aborting in file {} at line {}", file, line);
            mpi::abort(reason);
        } else {
            log::error_("Aborting in file {} at line {}", file, line);
            mpi::abort_(reason);
        }
    }

    template<typename lhs_t, typename rhs_t>
    static void compare_abort(const char *kind, const char *op, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym,
                              const char *rhs_sym,
                              const char *file, int line, bool collective, bool outcome, const std::string &reason) {
        if (collective || mpi::nrank() == 1) {
            outcome = mpi::all_land(outcome);
            if (!outcome) {
                log::error("{} {} failed", kind, op);
                log::error("LHS \"{}\" value is: {}", lhs_sym, convert::to_string(lhs));
                log::error("RHS \"{}\" value is: {}", rhs_sym, convert::to_string(rhs));
                abort(file, line, collective, reason);
            }
        } else {
            if (!outcome) {
                log::error_("{} {} failed", kind, op);
                log::error_("{} value is: {}", lhs_sym, convert::to_string(lhs));
                log::error_("{} value is: {}", rhs_sym, convert::to_string(rhs));
                abort(file, line, collective, reason);
            }
        }
    }

    template<typename T>
    static void compare_nearly_abort(const char *kind, const char *op, const T &lhs, const T &rhs,
                                     const arith::comp_t<T>& rtol, const arith::comp_t<T>& atol,
                                     const char *lhs_sym, const char *rhs_sym,
                                     const char *file, int line, bool collective, bool outcome, const std::string &reason) {
        if (collective || mpi::nrank() == 1) {
            outcome = mpi::all_land(outcome);
            if (!outcome) {
                log::error("{} {} failed with rtol={}, atol={}", kind, op, rtol, atol);
                log::error("LHS \"{}\" value is: {}", lhs_sym, convert::to_string(lhs));
                log::error("RHS \"{}\" value is: {}", rhs_sym, convert::to_string(rhs));
                abort(file, line, collective, reason);
            }
        } else {
            if (!outcome) {
                log::error_("{} {} failed with rtol={}, atol={}", kind, op, rtol, atol);
                log::error_("LHS \"{}\" value is: {}", lhs_sym, convert::to_string(lhs));
                log::error_("RHS \"{}\" value is: {}", rhs_sym, convert::to_string(rhs));
                abort(file, line, collective, reason);
            }
        }
    }

    template<typename T>
    static void compare_zero_abort(const char *kind, const char *sym, const char *file, int line, bool collective,
                                   T v, arith::comp_t<T> atol, const std::string &reason) {
        auto outcome = fptol::nearly_zero(v, atol);
        if (collective || mpi::nrank() == 1) {
            outcome = mpi::all_land(outcome);
            if (!outcome) {
                log::error("{} zero failed with atol={}", kind, atol);
                log::error("\"{}\" value is: {}", sym, convert::to_string(v));
                abort(file, line, collective, reason);
            }
        } else {
            if (!outcome) {
                log::error_("{} zero failed with atol={}", kind, atol);
                log::error_("\"{}\" value is: {}", sym, convert::to_string(v));
                abort(file, line, collective, reason);
            }
        }
    }

    static void bool_abort(const char *kind, const char *sym, const char *file, int line, bool collective,
                           bool outcome, bool truth, const std::string &reason) {
        auto right = truth ? "true" : "false";
        auto wrong = truth ? "false" : "true";
        if (collective || mpi::nrank() == 1) {
            if (!truth) outcome = !outcome;
            outcome = mpi::all_land(outcome); // true only if all input outcomes across all ranks are false

            if (!outcome) {
                log::error("{} {} failed", kind, right);
                log::error("{} value is {}", sym, wrong);
                abort(file, line, collective, reason);
            }
        } else {
            if (!truth) outcome = !outcome;
            if (!outcome) {
                log::error_("{} {} failed", kind, right);
                log::error_("{} value is {}", sym, wrong);
                abort(file, line, collective, reason);
            }
        }
    }

    template<typename lhs_t, typename rhs_t>
    static void eq(const char *kind, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym, const char *rhs_sym,
                   const char *file, int line, bool collective, const std::string &reason) {
        compare_abort(kind, "equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs == rhs, reason);
    }

    template<typename lhs_t, typename rhs_t>
    static void ne(const char *kind, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym, const char *rhs_sym,
                   const char *file, int line, bool collective, const std::string &reason) {
        compare_abort(kind, "not equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs != rhs, reason);
    }

    template<typename T>
    static void nearly_eq(const char *kind, const T &lhs, const T &rhs,
                          const arith::comp_t<T> &rtol, const arith::comp_t<T> &atol,
                          const char *lhs_sym, const char *rhs_sym,
                          const char *file, int line, bool collective, const std::string &reason) {
        compare_nearly_abort(kind, "equal", lhs, rhs, rtol, atol,lhs_sym, rhs_sym, file, line,
                             collective, fptol::nearly_equal(lhs, rhs, rtol, atol), reason);
    }

    template<typename T>
    static void nearly_eq(const char *kind, const T &lhs, const T &rhs,
                          const char *lhs_sym, const char *rhs_sym,
                          const char *file, int line, bool collective, const std::string &reason) {
        nearly_eq(kind, lhs, rhs, fptol::default_rtol_near<T>(), fptol::default_atol_near<T>(),
                  lhs_sym, rhs_sym, file, line, collective, reason);
    }


    template<typename T>
    static void numeric_eq(const char *kind, const T &lhs, const T &rhs,
                           const char *lhs_sym, const char *rhs_sym,
                           const char *file, int line, bool collective, const std::string &reason) {
        nearly_eq(kind, lhs, rhs, 0.0, fptol::default_atol_num<T>(), lhs_sym, rhs_sym, file, line, collective, reason);
    }

    template<typename T>
    static void nearly_zero(const char *kind, T v, const arith::comp_t<T> &atol, const char *sym,
                            const char *file, int line, bool collective, const std::string &reason) {
        compare_zero_abort(kind, sym, file, line, collective, v, atol, true, reason);
    }

    template<typename T>
    static void numeric_zero(const char *kind, T v, const char *sym,
                            const char *file, int line, bool collective, const std::string &reason) {
        nearly_zero(kind, v, fptol::default_atol_num<T>(), sym, file, line, collective, v, atol, true, reason);
    }


    template<typename lhs_t, typename rhs_t>
    static void lt(const char *kind, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym, const char *rhs_sym,
                   const char *file, int line, bool collective, const std::string &reason) {
        compare_abort(kind, "less than", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs < rhs, reason);
    }

    template<typename lhs_t, typename rhs_t>
    static void le(const char *kind, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym, const char *rhs_sym,
                   const char *file, int line, bool collective, const std::string &reason) {
        compare_abort(kind, "less than or equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs <= rhs,
                      reason);
    }

    template<typename lhs_t, typename rhs_t>
    static void gt(const char *kind, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym, const char *rhs_sym,
                   const char *file, int line, bool collective, const std::string &reason) {
        compare_abort(kind, "greater than", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs > rhs, reason);
    }

    template<typename lhs_t, typename rhs_t>
    static void ge(const char *kind, const lhs_t &lhs, const rhs_t &rhs, const char *lhs_sym, const char *rhs_sym,
                   const char *file, int line, bool collective, const std::string &reason) {
        compare_abort(kind, "greater than or equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs >= rhs,
                      reason);
    }

    static void is_true(const char *kind, bool v, const char *sym, const char *file, int line,
                        bool collective, const std::string &reason) {
        bool_abort(kind, sym, file, line, collective, v, true, reason);
    }

    static void is_false(const char *kind, bool v, const char *sym, const char *file, int line,
                         bool collective, const std::string &reason) {
        bool_abort(kind, sym, file, line, collective, v, false, reason);
    }
}

/**
 * calling this macro on var suppresses unused variable warnings and expresses the intention to use var only in the
 * debug builds
 */
#ifndef NDEBUG
#define DEBUG_ONLY(var)
#else
#define DEBUG_ONLY(var) dtype::unused(var);
#endif

#define MPI_EQ_BASE(kind, lhs, rhs, collective, reason) asserts::eq(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_NE_BASE(kind, lhs, rhs, collective, reason) asserts::ne(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);

#define MPI_NEARLY_EQ_TOL_BASE(kind, lhs, rhs, rtol, atol, collective, reason) asserts::nearly_eq(kind, lhs, rhs, rtol, atol, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_NEARLY_EQ_BASE(kind, lhs, rhs, collective, reason) asserts::nearly_eq(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_NUM_EQ_BASE(kind, lhs, rhs, collective, reason) asserts::numeric_eq(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);

#define MPI_NEARLY_ZERO_TOL_BASE(kind, v, collective, reason) asserts::nearly_zero(kind, v, atol, #v, __FILE__, __LINE__, collective, reason);
#define MPI_NEARLY_ZERO_BASE(kind, v, collective, reason) asserts::nearly_zero(kind, v, #v, __FILE__, __LINE__, collective, reason);
#define MPI_NUM_ZERO_BASE(kind, v, collective, reason) asserts::numeric_zero(kind, v, #v, __FILE__, __LINE__, collective, reason);

#define MPI_LT_BASE(kind, lhs, rhs, collective, reason) asserts::lt(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_LE_BASE(kind, lhs, rhs, collective, reason) asserts::le(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_GT_BASE(kind, lhs, rhs, collective, reason) asserts::gt(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_GE_BASE(kind, lhs, rhs, collective, reason) asserts::ge(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective, reason);
#define MPI_TRUE_BASE(kind, v, collective, reason) asserts::is_true(kind, v, #v, __FILE__, __LINE__, collective, reason);
#define MPI_FALSE_BASE(kind, v, collective, reason) asserts::is_false(kind, v, #v, __FILE__, __LINE__, collective, reason);
#define MPI_ABORT_BASE(collective, reason) asserts::abort(__FILE__, __LINE__, collective, reason);

#ifndef NDEBUG
#define DEBUG_ASSERT_EQ_ALL(lhs, rhs, reason) MPI_EQ_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_EQ(lhs, rhs, reason) MPI_EQ_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_NE_ALL(lhs, rhs, reason) MPI_NE_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_NE(lhs, rhs, reason) MPI_NE_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_NEARLY_EQ_TOL_ALL(lhs, rhs, rtol, atol, reason) MPI_NEARLY_EQ_TOL_BASE("ASSERT", lhs, rhs, rtol, atol, true, reason)
#define DEBUG_ASSERT_NEARLY_EQ_ALL(lhs, rhs, reason) MPI_NEARLY_EQ_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_NUM_EQ_ALL(lhs, rhs, reason) MPI_NUM_EQ_BASE("ASSERT", lhs, rhs, true, reason)

#define DEBUG_ASSERT_NEARLY_EQ_TOL(lhs, rhs, rtol, atol, reason) MPI_NEARLY_EQ_TOL_BASE("ASSERT", lhs, rhs, rtol, atol, false, reason)
#define DEBUG_ASSERT_NEARLY_EQ(lhs, rhs, reason) MPI_NEARLY_EQ_BASE("ASSERT", lhs, rhs, false, reason)
#define DEBUG_ASSERT_NUM_EQ(lhs, rhs, reason) MPI_NUM_EQ_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_NEARLY_ZERO_TOL_ALL(v, atol, reason) MPI_NEARLY_ZERO_TOL_BASE("ASSERT", v, atol, true, reason)
#define DEBUG_ASSERT_NEARLY_ZERO_ALL(v, atol, reason) MPI_NEARLY_ZERO_BASE("ASSERT", v, atol, true, reason)
#define DEBUG_ASSERT_NUM_ZERO_ALL(v, atol, reason) MPI_NUM_ZERO_BASE("ASSERT", v, true, reason)

#define DEBUG_ASSERT_NEARLY_ZERO_TOT(v, atol, reason) MPI_NEARLY_ZERO_TOL_BASE("ASSERT", v, atol, false, reason)
#define DEBUG_ASSERT_NEARLY_ZERO(v, atol, reason) MPI_NEARLY_ZERO_BASE("ASSERT", v, atol, false, reason)
#define DEBUG_ASSERT_NUM_ZERO(v, atol, reason) MPI_NUM_ZERO_BASE("ASSERT", v, false, reason)

#define DEBUG_ASSERT_NE_ALL(lhs, rhs, reason) MPI_NE_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_NE(lhs, rhs, reason) MPI_NE_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_LT_ALL(lhs, rhs, reason) MPI_LT_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_LT(lhs, rhs, reason) MPI_LT_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_LE_ALL(lhs, rhs, reason) MPI_LE_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_LE(lhs, rhs, reason) MPI_LE_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_GT_ALL(lhs, rhs, reason) MPI_GT_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_GT(lhs, rhs, reason) MPI_GT_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_GE_ALL(lhs, rhs, reason) MPI_GE_BASE("ASSERT", lhs, rhs, true, reason)
#define DEBUG_ASSERT_GE(lhs, rhs, reason) MPI_GE_BASE("ASSERT", lhs, rhs, false, reason)

#define DEBUG_ASSERT_TRUE_ALL(v, reason) MPI_TRUE_BASE("ASSERT", v, true, reason)
#define DEBUG_ASSERT_TRUE(v, reason) MPI_TRUE_BASE("ASSERT", v, false, reason)

#define DEBUG_ASSERT_FALSE_ALL(v, reason) MPI_FALSE_BASE("ASSERT", v, true, reason)
#define DEBUG_ASSERT_FALSE(v, reason) MPI_FALSE_BASE("ASSERT", v, false, reason)
#else
#define DEBUG_ASSERT_EQ_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_EQ(lhs, rhs, reason)

#define DEBUG_ASSERT_NE_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_NE(lhs, rhs, reason)

#define DEBUG_ASSERT_NEARLY_EQ_TOL_ALL(lhs, rhs, rtol, atol, reason)
#define DEBUG_ASSERT_NEARLY_EQ_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_NUM_EQ_ALL(lhs, rhs, reason)

#define DEBUG_ASSERT_NEARLY_EQ_TOL(lhs, rhs, rtol, atol, reason)
#define DEBUG_ASSERT_NEARLY_EQ(lhs, rhs, reason)
#define DEBUG_ASSERT_NUM_EQ(lhs, rhs, reason)

#define DEBUG_ASSERT_NEARLY_ZERO_TOL_ALL(v, atol, reason)
#define DEBUG_ASSERT_NEARLY_ZERO_ALL(v, atol, reason)
#define DEBUG_ASSERT_NUM_ZERO_ALL(v, atol, reason)

#define DEBUG_ASSERT_NEARLY_ZERO_TOT(v, atol, reason)
#define DEBUG_ASSERT_NEARLY_ZERO(v, atol, reason)
#define DEBUG_ASSERT_NUM_ZERO(v, atol, reason)

#define DEBUG_ASSERT_LT_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_LT(lhs, rhs, reason)

#define DEBUG_ASSERT_LE_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_LE(lhs, rhs, reason)

#define DEBUG_ASSERT_GT_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_GT(lhs, rhs, reason)

#define DEBUG_ASSERT_GE_ALL(lhs, rhs, reason)
#define DEBUG_ASSERT_GE(lhs, rhs, reason)

#define DEBUG_ASSERT_TRUE_ALL(v, reason)
#define DEBUG_ASSERT_TRUE(v, reason)

#define DEBUG_ASSERT_FALSE_ALL(v, reason)
#define DEBUG_ASSERT_FALSE(v, reason)
#endif

#define REQUIRE_EQ_ALL(lhs, rhs, reason) MPI_EQ_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_EQ(lhs, rhs, reason) MPI_EQ_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_NE_ALL(lhs, rhs, reason) MPI_NE_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_NE(lhs, rhs, reason) MPI_NE_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_NEARLY_EQ_TOL_ALL(lhs, rhs, rtol, atol, reason) MPI_NEARLY_EQ_TOL_BASE("REQUIRE", lhs, rhs, rtol, atol, true, reason)
#define REQUIRE_NEARLY_EQ_ALL(lhs, rhs, reason) MPI_NEARLY_EQ_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_NUM_EQ_ALL(lhs, rhs, reason) MPI_NUM_EQ_BASE("REQUIRE", lhs, rhs, true, reason)

#define REQUIRE_NEARLY_EQ_TOL(lhs, rhs, rtol, atol, reason) MPI_NEARLY_EQ_TOL_BASE("REQUIRE", lhs, rhs, rtol, atol, false, reason)
#define REQUIRE_NEARLY_EQ(lhs, rhs, reason) MPI_NEARLY_EQ_BASE("REQUIRE", lhs, rhs, false, reason)
#define REQUIRE_NUM_EQ(lhs, rhs, reason) MPI_NUM_EQ_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_NEARLY_ZERO_TOL_ALL(v, atol, reason) MPI_NEARLY_ZERO_TOL_BASE("REQUIRE", v, atol, true, reason)
#define REQUIRE_NEARLY_ZERO_ALL(v, atol, reason) MPI_NEARLY_ZERO_BASE("REQUIRE", v, atol, true, reason)
#define REQUIRE_NUM_ZERO_ALL(v, atol, reason) MPI_NUM_ZERO_BASE("REQUIRE", v, true, reason)

#define REQUIRE_NEARLY_ZERO_TOT(v, atol, reason) MPI_NEARLY_ZERO_TOL_BASE("REQUIRE", v, atol, false, reason)
#define REQUIRE_NEARLY_ZERO(v, atol, reason) MPI_NEARLY_ZERO_BASE("REQUIRE", v, atol, false, reason)
#define REQUIRE_NUM_ZERO(v, atol, reason) MPI_NUM_ZERO_BASE("REQUIRE", v, false, reason)

#define REQUIRE_LT_ALL(lhs, rhs, reason) MPI_LT_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_LT(lhs, rhs, reason) MPI_LT_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_LE_ALL(lhs, rhs, reason) MPI_LE_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_LE(lhs, rhs, reason) MPI_LE_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_GT_ALL(lhs, rhs, reason) MPI_GT_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_GT(lhs, rhs, reason) MPI_GT_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_GE_ALL(lhs, rhs, reason) MPI_GE_BASE("REQUIRE", lhs, rhs, true, reason)
#define REQUIRE_GE(lhs, rhs, reason) MPI_GE_BASE("REQUIRE", lhs, rhs, false, reason)

#define REQUIRE_TRUE_ALL(v, reason) MPI_TRUE_BASE("REQUIRE", v, true, reason)
#define REQUIRE_TRUE(v, reason) MPI_TRUE_BASE("REQUIRE", v, false, reason)

#define REQUIRE_FALSE_ALL(v, reason) MPI_FALSE_BASE("REQUIRE", v, true, reason)
#define REQUIRE_FALSE(v, reason) MPI_FALSE_BASE("REQUIRE", v, false, reason)

#define ABORT_ALL(reason) MPI_ABORT_BASE(true, reason)
#define ABORT(reason) MPI_ABORT_BASE(false, reason)

#endif //M7_MPIASSERT_H
