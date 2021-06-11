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
    ddx_sse2(
#elif defined(ISPC_TARGET_SSE4)
    void
    ddx_sse4(
#elif defined(ISPC_TARGET_AVX)
    void
    ddx_avx(
#elif defined(ISPC_TARGET_AVX2)
    void
    ddx_avx2(
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
    foreach (iz = 0 ... Nz, iy = 0 ... Ny, ix = 0 ... Nx) {
        dst[ix + x_pitch_dst * iy + xy_pitch_dst * iz] =
            fac * (src[(ix - 1) + x_pitch_src * iy + xy_pitch_src * iz] +
                   src[(ix + 1) + x_pitch_src * iy + xy_pitch_src * iz] +
                   src[ix + x_pitch_src * (iy - 1) + xy_pitch_src * iz] +
                   src[ix + x_pitch_src * (iy + 1) + xy_pitch_src * iz] +
                   src[ix + x_pitch_src * iy + xy_pitch_src * (iz - 1)] +
                   src[ix + x_pitch_src * iy + xy_pitch_src * (iz + 1)] -
                   cfac * src[ix + x_pitch_src * iy + xy_pitch_src * iz]);
    }
}
