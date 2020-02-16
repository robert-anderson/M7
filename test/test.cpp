//
// Created by Robert John Anderson on 2020-01-04.
//

#include <gtest/gtest.h>
#include "src/parallel/MPIWrapper.h"
#include "gtest/gtest.h"

int main(int argc, char **argv) {

    int out = 0;

    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    mpi::initialize(&argc, &argv);

    ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();
    if (!mpi::i_am_root()) {
        delete listeners.Release(listeners.default_result_printer());
    }

#if USE_MPI

#endif
    // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
    // pass and 1 if some test fails.
    auto result = RUN_ALL_TESTS();
    if (mpi::i_am_root()) out=result;

    mpi::finalize();
    return out;

}
