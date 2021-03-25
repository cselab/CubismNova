// File       : Laplacian.cpp
// Created    : Wed Mar 24 2021 04:30:24 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Compute the Laplacian on a Cartesian grid
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Math.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/IO/CartesianHDF.h"
#include "../../../Utils/Timer.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace Cubism;
using Utils::Timer;

// #define _EXPLICIT_LAB_LOAD_

int main(int argc, char *argv[])
{
    // Convenience data types
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using Point = typename Mesh::PointType;
    using MIndex = typename Mesh::MultiIndex;
    using ScalarGrid = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>;
    using DataType = typename ScalarGrid::DataType;
    using ScalarFieldType = typename ScalarGrid::BaseType;
    using FieldLab = Block::FieldLab<ScalarFieldType>;
    using Stencil = typename FieldLab::StencilType;

    // grid blocks and cells per block
    const MIndex nblocks((2 == argc) ? std::atoi(argv[1]) : 8);
    const MIndex block_cells(32);
    Timer t;
    ScalarGrid u(nblocks, block_cells);
    ScalarGrid lap_u(nblocks, block_cells); // memory to write back
    double t0 = t.stop();
    const Point h = u.getMesh().getCellSize(0); // uniform cells

    // initial condition
    auto u0 = [](ScalarFieldType &field) {
        const Mesh &bm =
            *field.getState().mesh; // block field mesh (partial mesh)
        for (auto &ci : bm[EntityType::Cell]) {
            const Point x = bm.getCoordsCell(ci); // cell center coordinate
            field[ci] = std::sin(2.0 * M_PI * x[0]) *
                        std::cos(6.0 * M_PI * x[1]) *
                        std::sin(4.0 * M_PI * x[2]);
        }
    };

    // initialize grid: iterate over all block fields contained in grid
    t.start();
    for (auto bf : u) {
        u0(*bf); // bf is an iterator, need to dereference it
    }
    t0 += t.stop();

    // dump initial grid into HDF5 file (note that data in `u` is double
    // precision, we dump the file in single precision)
    t.start();
    IO::CartesianWriteHDF<float>("init", "U", u, 0);
    double td = t.stop();

    // setup lab
    t.start();
    FieldLab lab;
    const Stencil s(-1, 2, false); // stencil (non-tensorial)
    // one lab per thread (use only one here)
    lab.allocate(s, u[0].getIndexRange());

    const MIndex ix{1, 0, 0};
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    const DataType ihx2 = 1. / (h[0] * h[0]);
    const DataType ihy2 = 1. / (h[1] * h[1]);
    const DataType ihz2 = 1. / (h[2] * h[2]);
#ifdef _EXPLICIT_LAB_LOAD_
    // explicit block index functor (getting it from ScalarGrid results in
    // periodic indexer)
    auto findex = u.getIndexFunctor();

    // process block fields in grid
    for (auto f : u) {
        // field reference for which we load the laboratory
        const ScalarFieldType &bf = *f;
        const MIndex &bi = bf.getState().block_index; // block index in grid

        // load the data into the previously allocated lab explicitly
        lab.loadData(bi, findex);

        // obtain block field where we write the result
        auto &res = lap_u[bi];

        // process lab (compute Laplacian)
        for (auto &i : bf.getIndexRange()) {
            // clang-format off
            const DataType ddx = ihx2 * (lab[i + ix] - 2 * lab[i] + lab[i - ix]);
            const DataType ddy = ihy2 * (lab[i + iy] - 2 * lab[i] + lab[i - iy]);
            const DataType ddz = ihz2 * (lab[i + iz] - 2 * lab[i] + lab[i - iz]);
            // clang-format on
            res[i] = ddx + ddy + ddz;
        }
    }
#else
    // process block fields in grid
    for (auto f : u) {
        // field reference for which we load the laboratory
        const ScalarFieldType &bf = *f;
        const MIndex &bi = bf.getState().block_index; // block index in grid

        // load the data into the previously allocated lab via the grid type
        // (performs additional checks in debug mode, preferred method)
        u.loadLab(bf, lab);

        // obtain block field where we write the result
        auto &res = lap_u[bi];

        // process lab (compute Laplacian)
        for (auto &i : bf.getIndexRange()) {
            // clang-format off
            const DataType ddx = ihx2 * (lab[i + ix] - 2 * lab[i] + lab[i - ix]);
            const DataType ddy = ihy2 * (lab[i + iy] - 2 * lab[i] + lab[i - iy]);
            const DataType ddz = ihz2 * (lab[i + iz] - 2 * lab[i] + lab[i - iz]);
            // clang-format on
            res[i] = ddx + ddy + ddz;
        }
    }
#endif /* _EXPLICIT_LAB_LOAD_ */
    const double t1 = t.stop();

    t.start();
    IO::CartesianWriteHDF<float>("laplacian", "div(grad(U))", lap_u, 0);
    td += t.stop();
    printf("ncells:  %d\ninit:    %e\ncompute: %e\ndump:    %e\ntotal:   %e\n",
           static_cast<int>(nblocks.prod()),
           t0,
           t1,
           td,
           t0 + t1 + td);

    return 0;
}
