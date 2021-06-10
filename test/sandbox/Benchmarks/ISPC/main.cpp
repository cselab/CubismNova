// File       : main.cpp
// Created    : Wed Jun 09 2021 10:51:40 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Main ISPC driver code
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Common.h"

#include "CFD_2ndOrder/Benchmark.h"

#include <cstdio>

// number of samples to collect
#ifndef _NSAMPLES_
#define _NSAMPLES_ 50
#endif /* _NSAMPLES_ */
constexpr int n_samples = _NSAMPLES_;

int main(int argc, char *argv[])
{
#ifdef _SINGLE_PRECISION_
    const char prec[] = "Single";
#else
    const char prec[] = "Double";
#endif /* _SINGLE_PRECISION_ */
    printf("BENCHMARK PRECISION: %s\n", prec);

    const int elements_per_dim = (2 == argc) ? std::atoi(argv[1]) : 128;

    // Benchmark objects
    CFD::Order2::Benchmark benchmark(n_samples, elements_per_dim);
    benchmark.run();

    return 0;
}
