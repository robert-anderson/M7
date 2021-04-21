//
// Created by rja on 29/01/2021.
//

#ifndef M7_MPIASSERT_H
#define M7_MPIASSERT_H

#include "MPIWrapper.h"
#include "src/core/io/Logging.h"

/**
 * these macros enable output of invoking file and line before quitting.
 * with a suffix "_", this is a forced quit on a subset of the MPI_COMM_WORLD
 * without this suffix, the quit is done collectively, which is preferred.
 * the former version should only be used inside selective functions which are
 * only invoked on a subset of MPI ranks
 */
#define MPI_ABORT_(msg) {log::error_("MPI_ABORT_ called at {}:{}", __FILE__, __LINE__); mpi::abort_(msg);}
#define MPI_ABORT(msg) {log::error("MPI_ABORT called at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}

//#define MPI_STOP_ALL(msg) {log::error("MPI_STOP_ALL called at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}


namespace asserts {

    template<typename lhs_t, typename rhs_t>
    static void compare_abort(const char* kind, const char* op, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                              const char* file, int line, bool collective, bool outcome){
        if (collective){
            outcome = mpi::all_land(outcome);
            if (!outcome) {
                log::error("{} {} failed in file {} at line {}", kind, op, file, line);
                log::error("LHS \"{}\" value is: {}", lhs_sym, lhs);
                log::error("RHS \"{}\" value is: {}", rhs_sym, rhs);
                mpi::abort(log::format("{} {} failure", kind, op));
            }
        }
        else {
            if (!outcome) {
                log::error_("{} {} failed", kind, op);
                log::error_("{} value is: {}", lhs_sym, lhs);
                log::error_("{} value is: {}", rhs_sym, rhs);
                mpi::abort_(log::format("{} {} failure", kind, op));
            }
        }
    }

    static void truth_abort(const char* kind, const char* sym,
                              const char* file, int line, bool collective, bool outcome){
        if (collective){
            outcome = mpi::all_land(outcome);
            if (!outcome) {
                log::error("{} true failed in file {} at line {}", kind, file, line);
                log::error("{} value is false", sym);
                mpi::abort(log::format("{} failure", kind));
            }
        }
        else {
            if (!outcome) {
                log::error_("{} true failed in file {} at line {}", kind, file, line);
                log::error_("{} value is false", sym);
                mpi::abort_(log::format("{} failure", kind));
            }
        }
    }

    template<typename lhs_t, typename rhs_t>
    static void eq(const char* kind, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                   const char* file, int line, bool collective){
        compare_abort(kind, "equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs==rhs);
    }

    template<typename lhs_t, typename rhs_t>
    static void ne(const char* kind, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                   const char* file, int line, bool collective){
        compare_abort(kind, "not equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs!=rhs);
    }

    template<typename lhs_t, typename rhs_t>
    static void lt(const char* kind, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                   const char* file, int line, bool collective){
        compare_abort(kind, "less than", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs<rhs);
    }

    template<typename lhs_t, typename rhs_t>
    static void le(const char* kind, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                   const char* file, int line, bool collective){
        compare_abort(kind, "less than or equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs<=rhs);
    }

    template<typename lhs_t, typename rhs_t>
    static void gt(const char* kind, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                   const char* file, int line, bool collective){
        compare_abort(kind, "greater than", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs>rhs);
    }

    template<typename lhs_t, typename rhs_t>
    static void ge(const char* kind, const lhs_t& lhs, const rhs_t& rhs, const char* lhs_sym, const char* rhs_sym,
                   const char* file, int line, bool collective){
        compare_abort(kind, "greater than or equal", lhs, rhs, lhs_sym, rhs_sym, file, line, collective, lhs>=rhs);
    }

    static void truth(const char* kind, bool v, const char* sym,
                   const char* file, int line, bool collective){
        truth_abort(kind, sym, file, line, collective, v);
    }
}

// suppresses unused variable warnings
#define TOUCH(var) (void)(var)

#define MPI_EQ_BASE(kind, lhs, rhs, collective) asserts::eq(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective);
#define MPI_NE_BASE(kind, lhs, rhs, collective) asserts::ne(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective);
#define MPI_LT_BASE(kind, lhs, rhs, collective) asserts::lt(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective);
#define MPI_LE_BASE(kind, lhs, rhs, collective) asserts::le(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective);
#define MPI_GT_BASE(kind, lhs, rhs, collective) asserts::gt(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective);
#define MPI_GE_BASE(kind, lhs, rhs, collective) asserts::ge(kind, lhs, rhs, #lhs, #rhs, __FILE__, __LINE__, collective);
#define MPI_TRUE_BASE(kind, v, collective) asserts::truth(kind, v, #v, __FILE__, __LINE__, collective);

#ifndef DNDEBUG
#define MPI_ASSERT_EQ(lhs, rhs) MPI_EQ_BASE("ASSERT", lhs, rhs, true)
#define MPI_ASSERT_EQ_(lhs, rhs) MPI_EQ_BASE("ASSERT", lhs, rhs, false)

#define MPI_ASSERT_NE(lhs, rhs) MPI_NE_BASE("ASSERT", lhs, rhs, true)
#define MPI_ASSERT_NE_(lhs, rhs) MPI_NE_BASE("ASSERT", lhs, rhs, false)

#define MPI_ASSERT_LT(lhs, rhs) MPI_LT_BASE("ASSERT", lhs, rhs, true)
#define MPI_ASSERT_LT_(lhs, rhs) MPI_LT_BASE("ASSERT", lhs, rhs, false)

#define MPI_ASSERT_LE(lhs, rhs) MPI_LE_BASE("ASSERT", lhs, rhs, true)
#define MPI_ASSERT_LE_(lhs, rhs) MPI_LE_BASE("ASSERT", lhs, rhs, false)

#define MPI_ASSERT_GT(lhs, rhs) MPI_GT_BASE("ASSERT", lhs, rhs, true)
#define MPI_ASSERT_GT_(lhs, rhs) MPI_GT_BASE("ASSERT", lhs, rhs, false)

#define MPI_ASSERT_GE(lhs, rhs) MPI_GE_BASE("ASSERT", lhs, rhs, true)
#define MPI_ASSERT_GE_(lhs, rhs) MPI_GE_BASE("ASSERT", lhs, rhs, false)

#define MPI_ASSERT_TRUE(v) MPI_TRUE_BASE("ASSERT", v, true)
#define MPI_ASSERT_TRUE_(v) MPI_TRUE_BASE("ASSERT", v, false)
#else

#define MPI_ASSERT_EQ(lhs, rhs) TOUCH(lhs); TOUCH(rhs);
#define MPI_ASSERT_EQ_(lhs, rhs) TOUCH(lhs); TOUCH(rhs);

#define MPI_ASSERT_NE(lhs, rhs) TOUCH(lhs); TOUCH(rhs);
#define MPI_ASSERT_NE_(lhs, rhs) TOUCH(lhs); TOUCH(rhs);

#define MPI_ASSERT_LT(lhs, rhs) TOUCH(lhs); TOUCH(rhs);
#define MPI_ASSERT_LT_(lhs, rhs) TOUCH(lhs); TOUCH(rhs);

#define MPI_ASSERT_LE(lhs, rhs) TOUCH(lhs); TOUCH(rhs);
#define MPI_ASSERT_LE_(lhs, rhs) TOUCH(lhs); TOUCH(rhs);

#define MPI_ASSERT_GT(lhs, rhs) TOUCH(lhs); TOUCH(rhs);
#define MPI_ASSERT_GT_(lhs, rhs) TOUCH(lhs); TOUCH(rhs);

#define MPI_ASSERT_GE(lhs, rhs) TOUCH(lhs); TOUCH(rhs);
#define MPI_ASSERT_GE_(lhs, rhs) TOUCH(lhs); TOUCH(rhs);

#define MPI_ASSERT_TRUE(v) TOUCH(v);
#define MPI_ASSERT_TRUE_(v) TOUCH(v);
#endif


#ifndef DNDEBUG
#define MPI_ASSERT(e, msg) {if (!(e)){log::error("MPI_ASSERT failed at {}:{}", __FILE__, __LINE__); mpi::abort_(msg);}}
#define MPI_ASSERT_ALL(e, msg) {if (!mpi::all_land((bool)(e))){log::error("MPI_ASSERT_ALL failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}}
#define MPI_ASSERT_ANY(e, msg) {if (!mpi::all_lor((bool)(e))){log::error("MPI_ASSERT_ANY failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}}
//#define MPI_ASSERT_ONE(e, irank, msg) {if (mpi::i_am(irank) && !e){log::error("MPI_ASSERT_ONE with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg);}}
//#define MPI_ASSERT_ALL_OTHERS(e, irank, msg) {if (!mpi::all_land(!mpi::i_am(irank)&&e)){log::error("MPI_ASSERT_ALL_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg);}}
//#define MPI_ASSERT_ANY_OTHERS(e, irank, msg) {if (!mpi::all_lor(!mpi::i_am(irank)&&e)){log::error("MPI_ASSERT_ANY_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg);}}
#else
#define MPI_ASSERT_ALL(e, msg)
#define MPI_ASSERT_ANY(e, msg)
#define MPI_ASSERT_ONE(e, irank, msg)
#define MPI_ASSERT_ALL_OTHERS(e, irank, msg)
#define MPI_ASSERT_ANY_OTHERS(e, irank, msg)
#endif

// REQUIRE macros perform the same role as ASSERTs, but are also defined in release build
#define MPI_REQUIRE(e, msg) {if (!(e)){log::error("MPI_REQUIRE failed at {}:{}", __FILE__, __LINE__); mpi::abort_(msg);}}
#define MPI_REQUIRE_ALL(e, msg) {if (!mpi::all_land((bool)(e))){log::error("MPI_REQUIRE_ALL failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}}
#define MPI_REQUIRE_ANY(e, msg) {if (!mpi::all_lor((bool)(e))){log::error("MPI_REQUIRE_ANY failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg);}}
//#define MPI_REQUIRE_ONE(e, irank, msg) {if (mpi::i_am(irank) && !e){log::error("MPI_REQUIRE_ONE with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg);}}
//#define MPI_REQUIRE_ALL_OTHERS(e, irank, msg) {if (!mpi::all_land(!mpi::i_am(irank)&&e)){log::error("MPI_REQUIRE_ALL_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg);}}
//#define MPI_REQUIRE_ANY_OTHERS(e, irank, msg) {if (!mpi::all_lor(!mpi::i_am(irank)&&e)){log::error("MPI_REQUIRE_ANY_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg);}}
#endif //M7_MPIASSERT_H




#define MPI_REQUIRE_EQ(lhs, rhs) MPI_EQ_BASE("REQUIRE", lhs, rhs, true)
#define MPI_REQUIRE_EQ_(lhs, rhs) MPI_EQ_BASE("REQUIRE", lhs, rhs, false)

#define MPI_REQUIRE_NE(lhs, rhs) MPI_NE_BASE("REQUIRE", lhs, rhs, true)
#define MPI_REQUIRE_NE_(lhs, rhs) MPI_NE_BASE("REQUIRE", lhs, rhs, false)

#define MPI_REQUIRE_LT(lhs, rhs) MPI_LT_BASE("REQUIRE", lhs, rhs, true)
#define MPI_REQUIRE_LT_(lhs, rhs) MPI_LT_BASE("REQUIRE", lhs, rhs, false)

#define MPI_REQUIRE_LE(lhs, rhs) MPI_LE_BASE("REQUIRE", lhs, rhs, true)
#define MPI_REQUIRE_LE_(lhs, rhs) MPI_LE_BASE("REQUIRE", lhs, rhs, false)

#define MPI_REQUIRE_GT(lhs, rhs) MPI_GT_BASE("REQUIRE", lhs, rhs, true)
#define MPI_REQUIRE_GT_(lhs, rhs) MPI_GT_BASE("REQUIRE", lhs, rhs, false)

#define MPI_REQUIRE_GE(lhs, rhs) MPI_GE_BASE("REQUIRE", lhs, rhs, true)
#define MPI_REQUIRE_GE_(lhs, rhs) MPI_GE_BASE("REQUIRE", lhs, rhs, false)

#define MPI_REQUIRE_TRUE(v) MPI_TRUE_BASE("REQUIRE", v, true)
#define MPI_REQUIRE_TRUE_(v) MPI_TRUE_BASE("REQUIRE", v, false)