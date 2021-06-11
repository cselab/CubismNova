// File       : ddx.ispc
// Created    : Thu Jun 10 2021 01:33:43 PM (+0200)
// Author     : Fabian Wermelinger
// Description: ISPC implementation for second derivative (ddx)
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifdef _SINGLE_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif /* _SINGLE_PRECISION_ */

export
#if defined(ISPC_TARGET_SSE2)
    void
    ddxSlice_sse2(
#elif defined(ISPC_TARGET_SSE4)
    void
    ddxSlice_sse4(
#elif defined(ISPC_TARGET_AVX)
    void
    ddxSlice_avx(
#elif defined(ISPC_TARGET_AVX2)
    void
    ddxSlice_avx2(
#else
#error "Unknown ISPC target platform"
#endif
        const uniform int Nx,
        const uniform int Ny,
        const uniform int Nz,
        const uniform Real *uniform src,
        const uniform int x_pitch_src,
        const uniform int xy_pitch_src,
        uniform Real *uniform dst,
        const uniform int x_pitch_dst,
        const uniform int xy_pitch_dst,
        const uniform Real factor)
{
    const uniform Real fac =
        1.0 / (factor * factor); // square of inverse grid spacing
    const uniform Real cfac = 6.0;

    const uniform Real *uniform s_zm;
    const uniform Real *uniform s_z0 = src - xy_pitch_src;
    const uniform Real *uniform s_zp = src;
    for (uniform int iz = 0; iz < Nz; ++iz) {
        // advance slices
        s_zm = s_z0;
        s_z0 = s_zp;
        s_zp = src + (iz + 1) * xy_pitch_src;
        Real *s_dst = dst + iz * xy_pitch_dst;

        // process xy dimensions
        foreach (iy = 0 ... Ny, ix = 0 ... Nx) {
            s_dst[ix + x_pitch_dst * iy] = s_z0[(ix - 1) + x_pitch_src * iy] +
                                           s_z0[(ix + 1) + x_pitch_src * iy] +
                                           s_z0[ix + x_pitch_src * (iy - 1)] +
                                           s_z0[ix + x_pitch_src * (iy + 1)];
        }

        // process z-dimension, scaling and center node
        foreach (iy = 0 ... Ny, ix = 0 ... Nx) {
            s_dst[ix + x_pitch_dst * iy] =
                fac *
                (s_dst[ix + x_pitch_dst * iy] + s_zm[ix + x_pitch_src * iy] +
                 s_zp[ix + x_pitch_src * iy] -
                 cfac * s_z0[ix + x_pitch_src * iy]);
        }
    }
}
