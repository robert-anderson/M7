//
// Created by Robert John Anderson on 2020-01-04.
//

#include <gtest/gtest.h> // googletest header file
#include <external/gtest-mpi-listener/include/gtest-mpi-listener.hpp>
#include "../src/parallel/MPIWrapper.h"

int main(int argc, char **argv) {

    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    MPIWrapper::initialize(&argc, &argv);

#if USE_MPI
    // Add object that will finalize MPI on exit; Google Test owns this pointer
    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);

    // Get the event listener list.
    ::testing::TestEventListeners &listeners =
            ::testing::UnitTest::GetInstance()->listeners();

    // Remove default listener: the default printer and the default XML printer
    ::testing::TestEventListener *l =
            listeners.Release(listeners.default_result_printer());

    // Adds MPI listener; Google Test owns this pointer
    listeners.Append(
            new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD)
    );
#endif
    // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
    // pass and 1 if some test fails.
    int result = RUN_ALL_TESTS();
    return 0;
}
