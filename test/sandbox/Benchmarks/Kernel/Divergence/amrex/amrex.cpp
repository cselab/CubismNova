// File       : amrex.cpp
// Created    : Tue Mar 03 2020 07:35:31 PM (-0800)
// Author     : Fabian Wermelinger
// Description: Compute divergence test kernel
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../../../Utils/Timer.h"

using Utils::Timer;

using namespace amrex;

void main_main(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    int one = 1;
    amrex::Initialize(one, argv);

    main_main(argc, argv);

    amrex::Finalize();
    return 0;
}

using Fab = FArrayBox;
using Vec = Vector<Real>;

template <typename TFab>
void finit(
    const Box &b, TFab &f, const Real *dx, const Real *start, const Real *end)
{
    const Real fac = 2 * M_PI;
    for (int iz = b.smallEnd(2); iz <= b.bigEnd(2); ++iz) {
        const Real z = start[2] + (0.5 + iz) * dx[2];
        for (int iy = b.smallEnd(1); iy <= b.bigEnd(1); ++iy) {
            const Real y = start[1] + (0.5 + iy) * dx[1];
            for (int ix = b.smallEnd(0); ix <= b.bigEnd(0); ++ix) {
                const Real x = start[0] + (0.5 + ix) * dx[0];
                f(ix, iy, iz, 0) =
                    std::sin(fac * x) * std::cos(fac * y) * std::sin(fac * z);
                f(ix, iy, iz, 1) =
                    std::cos(fac * x) * std::cos(fac * y) * std::sin(fac * z);
                f(ix, iy, iz, 2) =
                    std::cos(fac * x) * std::cos(fac * y) * std::sin(fac * z);
            }
        }
    }
}

void main_main(int argc, char *argv[])
{
    static_assert(AMREX_SPACEDIM == 3, "AMREX_SPACEDIM must be 3");
    // AMREX_SPACEDIM: number of dimensions
    const int max_grid_size = 32; // 32 cells per block
    Vector<int> is_periodic(AMREX_SPACEDIM,
                            1); // periodic in all direction by default
    const int nblocks = (2 == argc) ? std::atoi(argv[1]) : 8;
    const int n_cell = nblocks * max_grid_size;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
        IntVect dom_hi(AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size"
        // along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical box, [0,1] in each direction.
        RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                         {AMREX_D_DECL(1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());
    }

    // How Boxes are distributed among MPI processes
    Timer t;
    DistributionMapping dm(ba);
    MultiFab phi(ba, dm, 3, 1); // 3 components; 1 ghost cell
    MultiFab res(ba, dm, 1, 0); // 1 component; 0 ghosts
    double t0 = t.stop();

    // initial condition
    // MFIter = MultiFab Iterator
    t.start();
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto ary = phi.array(mfi);
        finit(bx, ary, geom.CellSize(), geom.ProbLo(), geom.ProbHi());
    }
    t0 += t.stop();

#ifdef _DUMP_
    t.start();
    {
        const std::string &pltfile = amrex::Concatenate("plt", 0, 5);
        WriteSingleLevelPlotfile(
            pltfile, phi, {"phi0", "phi1", "phi2"}, geom, 0.0, 0);
    }
    double td = t.stop();
#endif /* _DUMP_ */

    t.start();
    // Fill ghost cells
    phi.FillBoundary(geom.periodicity(), true);
    // process blocks
    const Real *dx = geom.CellSize();
    const Real ih0 = 1.0 / dx[0];
    const Real ih1 = 1.0 / dx[1];
    const Real ih2 = 1.0 / dx[2];
    for (MFIter mfi(res); mfi.isValid(); ++mfi) {
        const Box &b = mfi.validbox();
        const auto src = phi.array(mfi);
        const auto dst = res.array(mfi);
        for (int iz = b.smallEnd(2); iz <= b.bigEnd(2); ++iz) {
            for (int iy = b.smallEnd(1); iy <= b.bigEnd(1); ++iy) {
                for (int ix = b.smallEnd(0); ix <= b.bigEnd(0); ++ix) {
                    const Real dd0_dx =
                        src(ix + 1, iy, iz, 0) - src(ix - 1, iy, iz, 0);
                    const Real dd1_dy =
                        src(ix, iy + 1, iz, 1) - src(ix, iy - 1, iz, 1);
                    const Real dd2_dz =
                        src(ix, iy, iz + 1, 2) - src(ix, iy, iz - 1, 2);
                    dst(ix, iy, iz, 0) =
                        0.5 * (dd0_dx * ih0 + dd1_dy * ih1 + dd2_dz * ih2);
                }
            }
        }
    }
    const double t1 = t.stop();

#ifdef _DUMP_
    {
        t.start();
        const std::string &pltfile = amrex::Concatenate("plt", 1, 5);
        WriteSingleLevelPlotfile(pltfile, res, {"res"}, geom, 0.0, 0);
        td += t.stop();
    }
    printf("%d\t%e\t%e\t%e\t%e\n",
           nblocks * nblocks * nblocks,
           t0,
           t1,
           td,
           t0 + t1 + td);
#else
    printf("%d\t%e\t%e\t%e\t%e\n",
           nblocks * nblocks * nblocks,
           t0,
           t1,
           0.,
           t0 + t1);
#endif /* _DUMP_ */
}
