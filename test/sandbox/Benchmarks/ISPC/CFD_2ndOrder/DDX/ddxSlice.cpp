// File       : Gold.cpp
// Created    : Thu Jun 10 2021 09:57:19 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Gold reference implementations (2D slice processing)
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "ddxSlice.h"

inline void processSlice_xy(const int Nx,
                            const int Ny,
                            const Real *src,
                            const int src_pitch,
                            Real *dst,
                            const int dst_pitch)
{
    // 3 Flop
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            dst[ix + dst_pitch * iy] =
                src[ix - 1 + src_pitch * iy] + src[ix + 1 + src_pitch * iy] +
                src[ix + src_pitch * (iy - 1)] + src[ix + src_pitch * (iy + 1)];
        }
    }
}

inline void processSlice_z(const int Nx,
                           const int Ny,
                           const Real *zm,
                           const Real *z0,
                           const Real *zp,
                           const int src_pitch,
                           Real *dst,
                           const int dst_pitch,
                           const Real fac)
{
    // 5 Flop
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            dst[ix + dst_pitch * iy] =
                fac * (dst[ix + dst_pitch * iy] + zm[ix + src_pitch * iy] +
                       zp[ix + src_pitch * iy] -
                       static_cast<Real>(6) * z0[ix + src_pitch * iy]);
        }
    }
}

// second derivative
#if defined(__GNUC__) || defined(__clang__)
// no auto vectorization
__attribute__((optimize("no-tree-vectorize")))
#endif
void CFD::Order2::Gold::ddxSlice(const int Nx,
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

    // slice pointer
    const Real *s_zm = nullptr;
    const Real *s_z0 = src - xy_pitch_src;
    const Real *s_zp = src;
    for (int iz = 0; iz < Nz; ++iz) {
        // advance slices
        s_zm = s_z0;
        s_z0 = s_zp;
        s_zp = src + (iz + 1) * xy_pitch_src;
        Real *s_dst = dst + iz * xy_pitch_dst;

        // process xy dimensions
        processSlice_xy(Nx, Ny, s_z0, x_pitch_src, s_dst, x_pitch_dst);

        // process z-dimension, scaling and center node
        processSlice_z(
            Nx, Ny, s_zm, s_z0, s_zp, x_pitch_src, s_dst, x_pitch_dst, fac);
    }
}

// auto vectorized
void CFD::Order2::Gold::ddxSliceTreeVec(const int Nx,
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

    // slice pointer
    const Real *s_zm = nullptr;
    const Real *s_z0 = src - xy_pitch_src;
    const Real *s_zp = src;
    for (int iz = 0; iz < Nz; ++iz) {
        // advance slices
        s_zm = s_z0;
        s_z0 = s_zp;
        s_zp = src + (iz + 1) * xy_pitch_src;
        Real *s_dst = dst + iz * xy_pitch_dst;

        // process xy dimensions
        processSlice_xy(Nx, Ny, s_z0, x_pitch_src, s_dst, x_pitch_dst);

        // process z-dimension, scaling and center node
        processSlice_z(
            Nx, Ny, s_zm, s_z0, s_zp, x_pitch_src, s_dst, x_pitch_dst, fac);
    }
}
