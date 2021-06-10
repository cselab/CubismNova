// File       : Benchmark.h
// Created    : Thu Jun 10 2021 09:41:23 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Centered finite differences 2nd order
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef BENCHMARK_H_BBOGPIXL
#define BENCHMARK_H_BBOGPIXL

#include "BaseBenchmark.h"

namespace CFD
{
namespace Order2
{
class Benchmark : public BaseBenchmark
{
public:
    Benchmark(const int n_samples,
              const int n_elements_per_dim,
              const bool is_tensorial = false)
        // stencil start = -1, stencil end = 2 (exclusive)
        : BaseBenchmark(n_samples, n_elements_per_dim, -1, 2, is_tensorial)
    {
    }

    void run() override;
};
} // namespace Order2
} // namespace CFD

#endif /* BENCHMARK_H_BBOGPIXL */
