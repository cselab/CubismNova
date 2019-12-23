// File       : TestProfiler.cpp
// Created    : Mon Dec 23 2019 03:37:08 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Profiler sandbox tests
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Util/Profiler.h"

#include <chrono>
#include <mpi.h>
#include <thread>

using namespace Cubism::Util;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    Profiler p("MPI", MPI_COMM_WORLD);
    for (int i = 0; i < 5; ++i) {
        p.push("Outer Loop");
        for (int s = 0; s < 10; ++s) {
            p.push("Sleep 100ms");
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            p.pop();

            p.push("Sleep 200ms");
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            p.pop();
        }
        p.pop();
        p.printReport();
    }
    p.printReport();

    MPI_Finalize();
    return 0;
}
