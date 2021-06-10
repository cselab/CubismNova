// File       : Benchmark.h
// Created    : Thu Jun 10 2021 09:41:23 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Centered finite differences 2nd order
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "CFD_2ndOrder/Benchmark.h"
#include "CFD_2ndOrder/Kernels.h"
#include <cstdio>

namespace CFD
{
namespace Order2
{
void Benchmark::run()
{
    printf("RUNNING: CFD_2ndOrder\n");

    Result gold = this->benchmark_("GOLD_DDX", Gold::ddx, Gold::loop_flop_ddx);
    this->benchmark_(
        "TREEVEC_DDX", Gold::ddx_tree_vec, Gold::loop_flop_ddx_tree_vec, &gold);
}
} // namespace Order2
} // namespace CFD
