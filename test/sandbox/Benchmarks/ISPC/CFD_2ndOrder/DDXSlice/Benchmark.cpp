// File       : Benchmark.h
// Created    : Thu Jun 10 2021 09:41:23 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Centered finite differences 2nd order (sliced implementation)
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Benchmark.h"
#include "Kernels.h"

namespace CFD
{
namespace Order2
{
namespace DDXSlice
{
void Benchmark::run()
{
    printf("RUNNING: CFD_2ndOrder (2nd derivative, 2D slice processing)\n");

    Result gold =
        this->benchmark_("GOLD", Gold::ddxSlice, Gold::loop_flop_ddx_slice);

    this->benchmark_("AUTOVEC_GOLD",
                     Gold::ddxSliceTreeVec,
                     Gold::loop_flop_ddx_slice_tree_vec,
                     &gold);

    this->benchmark_(
        "ISPC_SSE2", ispc::ddxSlice_sse2, ispc::loop_flop_ddx_slice, &gold);

    this->benchmark_(
        "ISPC_SSE4", ispc::ddxSlice_sse4, ispc::loop_flop_ddx_slice, &gold);

    this->benchmark_(
        "ISPC_AVX", ispc::ddxSlice_avx, ispc::loop_flop_ddx_slice, &gold);

    this->benchmark_(
        "ISPC_AVX2", ispc::ddxSlice_avx2, ispc::loop_flop_ddx_slice, &gold);
}

} // namespace DDXSlice
} // namespace Order2
} // namespace CFD
