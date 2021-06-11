// File       : Kernels.h
// Created    : Thu Jun 10 2021 09:56:12 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Benchmark kernels
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef KERNELS_H_9BARLO6D
#define KERNELS_H_9BARLO6D

#include "Gold.h"

// ISPC kernels
namespace CFD
{
namespace Order2
{
#include "ddxSlice_ispc_avx.h"
#include "ddxSlice_ispc_avx2.h"
#include "ddxSlice_ispc_sse2.h"
#include "ddxSlice_ispc_sse4.h"
namespace ispc
{
constexpr int loop_flop_ddx_slice = 8;
} // namespace ispc
} // namespace Order2
} // namespace CFD

#endif /* KERNELS_H_9BARLO6D */
