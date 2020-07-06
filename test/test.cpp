//
// Created by Robert John Anderson on 2020-01-04.
//

#include <gtest/gtest.h>
#include <fstream>
#include "src/core/parallel/MPIWrapper.h"

int main(int argc, char **argv) {

    int out = 0;

    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    mpi::initialize(&argc, &argv);

    std::streambuf *original_stdout_buffer = nullptr;
    std::streambuf *original_stderr_buffer = nullptr;
    std::ofstream ofstdout, ofstderr;

    if (!mpi::i_am_root()) {
        /*
         * only allow standard output and error from the root MPI rank
         */
#ifdef DNDEBUG
        std::cout.setstate(std::ios_base::failbit);
        std::cerr.setstate(std::ios_base::failbit);
#else
        original_stdout_buffer = std::cout.rdbuf();
        ofstdout = std::ofstream("rank_"+std::to_string(mpi::irank())+".out");
        std::cout.rdbuf(ofstdout.rdbuf());

        original_stderr_buffer = std::cerr.rdbuf();
        ofstderr = std::ofstream("rank_"+std::to_string(mpi::irank())+".err");
        std::cerr.rdbuf(ofstderr.rdbuf());
#endif
    }

    ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();
    if (!mpi::i_am_root()) {
        delete listeners.Release(listeners.default_result_printer());
    }

    // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
    // pass and 1 if some test fails.
    auto result = RUN_ALL_TESTS();

    if (!mpi::i_am_root()) {
        /*
         * reinstate original buffers
         */
#ifndef DNDEBUG
        std::cout.rdbuf(original_stdout_buffer);
        std::cerr.rdbuf(original_stderr_buffer);
#endif
    }

    mpi::finalize();
    out = result;
    return out;
}
