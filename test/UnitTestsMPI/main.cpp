// File       : main.cpp
// Created    : Sat Jan 11 2020 11:41:24 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Main MPI test driver
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "gtest-mpi-listener.hpp"
#include "gtest/gtest.h"
#include <mpi.h>

int main(int argc, char **argv)
{
    // Filter out Google Test arguments
    ::testing::InitGoogleTest(&argc, argv);

    // Initialize MPI
    MPI_Init(&argc, &argv);

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
        new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));

    // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
    // pass and 1 if some test fails.
    return RUN_ALL_TESTS();
}
