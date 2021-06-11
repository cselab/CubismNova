// File       : Gold.cpp
// Created    : Thu Jun 10 2021 09:57:19 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Gold reference implementations
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Gold.h"

#define LOOP_BODY()                                                            \
    do {                                                                       \
        for (int iz = 0; iz < Nz; ++iz) {                                      \
            for (int iy = 0; iy < Ny; ++iy) {                                  \
                for (int ix = 0; ix < Nx; ++ix) {                              \
                    dst[ix + x_pitch_dst * iy + xy_pitch_dst * iz] =           \
                        fac *                                                  \
                        (src[(ix - 1) + x_pitch_src * iy +                     \
                             xy_pitch_src * iz] +                              \
                         src[(ix + 1) + x_pitch_src * iy +                     \
                             xy_pitch_src * iz] +                              \
                         src[ix + x_pitch_src * (iy - 1) +                     \
                             xy_pitch_src * iz] +                              \
                         src[ix + x_pitch_src * (iy + 1) +                     \
                             xy_pitch_src * iz] +                              \
                         src[ix + x_pitch_src * iy +                           \
                             xy_pitch_src * (iz - 1)] +                        \
                         src[ix + x_pitch_src * iy +                           \
                             xy_pitch_src * (iz + 1)] -                        \
                         static_cast<Real>(6) *                                \
                             src[ix + x_pitch_src * iy + xy_pitch_src * iz]);  \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

// second derivative
#if defined(__GNUC__) || defined(__clang__)
// no auto vectorization
__attribute__((optimize("no-tree-vectorize")))
#endif
void CFD::Order2::Gold::ddx(const int Nx,
                       const int Ny,
                       const int Nz,
                       const Real *src,
                       const int x_pitch_src,
                       const int xy_pitch_src,
                       Real *dst,
                       const int x_pitch_dst,
                       const int xy_pitch_dst,
                       const Real factor)
{
    const Real fac = 1.0 / (factor * factor); // square of inverse grid spacing
    LOOP_BODY();
}

// auto vectorized
void CFD::Order2::Gold::ddxTreeVec(const int Nx,
                                   const int Ny,
                                   const int Nz,
                                   const Real *src,
                                   const int x_pitch_src,
                                   const int xy_pitch_src,
                                   Real *dst,
                                   const int x_pitch_dst,
                                   const int xy_pitch_dst,
                                   const Real factor)
{
    const Real fac = 1.0 / (factor * factor); // square of inverse grid spacing
    LOOP_BODY();
}

#undef LOOP_BODY
