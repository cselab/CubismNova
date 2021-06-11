// File       : Gold.h
// Created    : Thu Jun 10 2021 09:57:19 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Gold reference (2D slice processing)
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef GOLD_H_3ZDE5UKL
#define GOLD_H_3ZDE5UKL

#include "Common.h"

namespace CFD
{
namespace Order2
{
namespace Gold
{

// second derivative
constexpr int loop_flop_ddx_slice = 8;
void ddxSlice(const int Nx,
              const int Ny,
              const int Nz,
              const Real *src,
              const int x_pitch_src,
              const int xy_pitch_src,
              Real *dst,
              const int x_pitch_dst,
              const int xy_pitch_dst,
              const Real factor);

constexpr int loop_flop_ddx_slice_tree_vec = 8;
void ddxSliceTreeVec(const int Nx,
                     const int Ny,
                     const int Nz,
                     const Real *src,
                     const int x_pitch_src,
                     const int xy_pitch_src,
                     Real *dst,
                     const int x_pitch_dst,
                     const int xy_pitch_dst,
                     const Real factor);
} // namespace Gold
} // namespace Order2
} // namespace CFD

#endif /* GOLD_H_3ZDE5UKL */
