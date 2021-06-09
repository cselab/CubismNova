// File       : main.cpp
// Created    : Wed Jun 09 2021 10:51:40 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Main ISPC driver code
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Benchmark.h"

// benchmark data type
using Real = float;

// number of samples to collect
constexpr int n_samples = 20;


int main(int argc, char *argv[])
{
    // Benchmark object
    const int elements_per_dim = (2 == argc) ? std::atoi(argv[1]) : 32;
    // FD::CenteredSecondOrder<Real> benchmark(n_samples, elements_per_dim);
    Benchmark<Real> benchmark(n_samples, elements_per_dim, -1, 2);

    // info float precision

    // benchmark.writeTestData();
    benchmark.run();

    return 0;
}
