//
// Created by rja on 29/01/2021.
//

#ifndef M7_MPIASSERT_H
#define M7_MPIASSERT_H

#include "MPIWrapper.h"
#include "src/core/io/Logging.h"

// these macros enable output of invoking file and line
#define MPI_ABORT_(msg) {log::error_("MPI_ABORT_ called at {}:{}", __FILE__, __LINE__); mpi::abort_(msg, 0.0);}
#define MPI_ABORT(msg) {log::error("MPI_ABORT called at {}:{}", __FILE__, __LINE__); mpi::abort(msg, 0.0);}
#define MPI_STOP_ALL(msg) {log::error("MPI_STOP_ALL called at {}:{}", __FILE__, __LINE__); mpi::abort(msg, 0.0);}


#ifndef DNDEBUG
#define MPI_ASSERT(e, msg) {if (!(e)){log::error("MPI_ASSERT failed at {}:{}", __FILE__, __LINE__); mpi::abort_(msg);}}
#define MPI_ASSERT_ALL(e, msg) {if (!mpi::all_land((bool)(e))){log::error("MPI_ASSERT_ALL failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg,0);}}
#define MPI_ASSERT_ANY(e, msg) {if (!mpi::all_lor((bool)(e))){log::error("MPI_ASSERT_ANY failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg,0);}}
//#define MPI_ASSERT_ONE(e, irank, msg) {if (mpi::i_am(irank) && !e){log::error("MPI_ASSERT_ONE with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg,0);}}
//#define MPI_ASSERT_ALL_OTHERS(e, irank, msg) {if (!mpi::all_land(!mpi::i_am(irank)&&e)){log::error("MPI_ASSERT_ALL_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg,0);}}
//#define MPI_ASSERT_ANY_OTHERS(e, irank, msg) {if (!mpi::all_lor(!mpi::i_am(irank)&&e)){log::error("MPI_ASSERT_ANY_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg,0);}}
#else
#define MPI_ASSERT_ALL(e, msg)
#define MPI_ASSERT_ANY(e, msg)
#define MPI_ASSERT_ONE(e, irank, msg)
#define MPI_ASSERT_ALL_OTHERS(e, irank, msg)
#define MPI_ASSERT_ANY_OTHERS(e, irank, msg)
#endif

// REQUIRE macros perform the same role as ASSERTs, but are also defined in release build
#define MPI_REQUIRE(e, msg) {if (!(e)){log::error("MPI_REQUIRE failed at {}:{}", __FILE__, __LINE__); mpi::abort_(msg);}}
#define MPI_REQUIRE_ALL(e, msg) {if (!mpi::all_land((bool)(e))){log::error("MPI_REQUIRE_ALL failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg,0);}}
#define MPI_REQUIRE_ANY(e, msg) {if (!mpi::all_lor((bool)(e))){log::error("MPI_REQUIRE_ANY failed at {}:{}", __FILE__, __LINE__); mpi::abort(msg,0);}}
//#define MPI_REQUIRE_ONE(e, irank, msg) {if (mpi::i_am(irank) && !e){log::error("MPI_REQUIRE_ONE with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg,0);}}
//#define MPI_REQUIRE_ALL_OTHERS(e, irank, msg) {if (!mpi::all_land(!mpi::i_am(irank)&&e)){log::error("MPI_REQUIRE_ALL_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg,0);}}
//#define MPI_REQUIRE_ANY_OTHERS(e, irank, msg) {if (!mpi::all_lor(!mpi::i_am(irank)&&e)){log::error("MPI_REQUIRE_ANY_OTHERS with irank={} failed at {}:{}", irank, __FILE__, __LINE__); mpi::abort(msg,0);}}
#endif //M7_MPIASSERT_H
