// File       : Benchmark.h
// Created    : Thu Jun 10 2021 09:41:23 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Centered finite differences 2nd order (sliced implementation)
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef BENCHMARK_H_PPTUF59K
#define BENCHMARK_H_PPTUF59K

#include "../DDX/Benchmark.h"

namespace CFD
{
namespace Order2
{
namespace DDXSlice
{
class Benchmark : public DDX::Benchmark
{
public:
    Benchmark(const int n_samples,
              const int n_elements_per_dim,
              const bool is_tensorial = false)
        : DDX::Benchmark(n_samples, n_elements_per_dim, is_tensorial)
    {
    }

    void run() override;
};
} // namespace DDX
} // namespace Order2
} // namespace CFD

#endif /* BENCHMARK_H_PPTUF59K */
