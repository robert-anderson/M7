//
// Created by Robert John Anderson on 2020-01-04.
//

#include <gtest/gtest.h>
#include <fstream>
#include <memory>
#include <M7_lib/io/Logging.h>
#include "M7_lib/parallel/MPIWrapper.h"

int main(int argc, char **argv) {
    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    mpi::initialize(&argc, &argv);
    log::initialize();

    //std::streambuf *original_stdout_buffer = nullptr;
    //std::streambuf *original_stderr_buffer = nullptr;
    std::unique_ptr<std::ofstream> ofstdout, ofstderr;
    if (!mpi::i_am_root()) {
        /*
         * only allow standard output and error from the root MPI rank
         */
        //std::cout.setstate(std::ios_base::failbit);
        //std::cerr.setstate(std::ios_base::failbit);
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
        //std::cout.rdbuf(original_stdout_buffer);
        //std::cerr.rdbuf(original_stderr_buffer);
    }

    mpi::finalize();
    return result;
}
