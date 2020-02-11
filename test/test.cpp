//
// Created by Robert John Anderson on 2020-01-04.
//

#include <gtest/gtest.h> // googletest header file
//#include <external/gtest-mpi-listener/include/gtest-mpi-listener.hpp>
#include "src/parallel/MPIWrapper.h"
#include "gtest/gtest.h"

int main(int argc, char **argv) {

    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    MPIWrapper::initialize(&argc, &argv);
    MPIWrapper mpi;

    ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();
    if (!mpi.i_am_root()) {
        delete listeners.Release(listeners.default_result_printer());
    }

#if USE_MPI
    /*
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
    );*/
#endif
    // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
    // pass and 1 if some test fails.
    auto result = RUN_ALL_TESTS();

    MPIWrapper::finalize();
    return result;

}
