// File       : test-cubismnova.cpp
// Created    : Mon Feb 17 2020 08:58:17 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Compute divergence test kernel
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/DataLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/IO/CartesianHDF.h"
#include "Cubism/Math.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include <cmath>

using namespace Cubism;

int main(void)
{
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using Point = typename Mesh::PointType;
    using MIndex = typename Mesh::MultiIndex;
    using VGrid = Grid::Cartesian<float, Mesh, EntityType::Cell, 1>;
    using SGrid = Grid::Cartesian<float, Mesh, EntityType::Cell, 0>;
    using DataType = typename VGrid::DataType;
    using FieldType = typename VGrid::BaseType;
    using Lab = Block::DataLab<typename FieldType::FieldType>;
    using Stencil = typename Lab::StencilType;
    // grid blocks and cells per block
    const MIndex nblocks(16);
    const MIndex block_cells(32);
    SGrid result(nblocks, block_cells);
    VGrid grid(nblocks, block_cells);
    const Point h = grid.getMesh().getCellSize(0);

    auto finit = [h](FieldType &b) {
        const DataType fac = 2 * M_PI;
        const Mesh &bm = *b.getState().mesh;
        // auto begin = bm[EntityType::Cell].begin();
        // auto end = bm[EntityType::Cell].end();
        // for (auto ci = begin; ci != end; ++ci) {
        //     const Point x = bm.getCoordsCell(*ci);
        //     b[0][*ci] = std::sin(fac * x[0]) * std::cos(fac * x[1]) *
        //                 std::sin(fac * x[2]);
        //     b[1][*ci] = std::cos(fac * x[0]) * std::cos(fac * x[1]) *
        //                 std::sin(fac * x[2]);
        //     b[2][*ci] = std::cos(fac * x[0]) * std::cos(fac * x[1]) *
        //                 std::sin(fac * x[2]);
        for (auto &ci : bm[EntityType::Cell]) {
            const Point x = bm.getCoordsCell(ci);
            // const Point x = bm.getCoords(ci, EntityType::Cell);
            b[0][ci] = std::sin(fac * x[0]) * std::cos(fac * x[1]) *
                       std::sin(fac * x[2]);
            b[1][ci] = std::cos(fac * x[0]) * std::cos(fac * x[1]) *
                       std::sin(fac * x[2]);
            b[2][ci] = std::cos(fac * x[0]) * std::cos(fac * x[1]) *
                       std::sin(fac * x[2]);
        }
    };

    // initialize grid
    for (auto f : grid) {
        finit(*f);
    }
#ifdef _DUMP_
    IO::CartesianWriteHDF<float>("init", "U", grid, 0);
#endif /* _DUMP_ */

    // setup lab
    Lab lab0, lab1, lab2;
    const Stencil s(-1, 2, false);     // stencil
    const MIndex sbegin(s.getBegin()); // stencil begin
    const MIndex send(s.getEnd() - 1); // stencil end
    lab0.allocate(s, grid[0].getIndexRange());
    lab1.allocate(s, grid[0].getIndexRange());
    lab2.allocate(s, grid[0].getIndexRange());

    // process block fields
    using IndexFunctor = typename VGrid::IndexFunctor;
    auto findex0 = grid.getIndexFunctor(0); // grid block index functor comp = 0
    auto findex1 = grid.getIndexFunctor(1); // grid block index functor comp = 1
    auto findex2 = grid.getIndexFunctor(2); // grid block index functor comp = 2

    const DataType ih0 = 1. / h[0];
    const DataType ih1 = 1. / h[1];
    const DataType ih2 = 1. / h[2];
    for (auto f : grid)
    {
        const FieldType &bf = *f; // reference for field
        lab0.loadData(bf.getState().block_index, findex0);
        lab1.loadData(bf.getState().block_index, findex1);
        lab2.loadData(bf.getState().block_index, findex2);

        auto &sf = result[bf.getState().block_index];

        // process lab
        const MIndex ix{1, 0, 0};
        const MIndex iy{0, 1, 0};
        const MIndex iz{0, 0, 1};
        // auto begin = bf.getIndexRange().begin();
        // auto end = bf.getIndexRange().end();
        // for (auto i = begin; i != end; ++i) {
        //     const DataType dd0_dx = lab0[*i + ix] - lab0[*i - ix];
        //     const DataType dd1_dy = lab1[*i + iy] - lab1[*i - iy];
        //     const DataType dd2_dz = lab2[*i + iz] - lab2[*i - iz];
        //     sf[*i] = 0.5 * (dd0_dx * ih0 + dd1_dy * ih1 + dd2_dz * ih2);
        for (auto i : bf.getIndexRange()) {
            const DataType dd0_dx = lab0[i + ix] - lab0[i - ix];
            const DataType dd1_dy = lab1[i + iy] - lab1[i - iy];
            const DataType dd2_dz = lab2[i + iz] - lab2[i - iz];
            sf[i] = 0.5 * (dd0_dx * ih0 + dd1_dy * ih1 + dd2_dz * ih2);
        }
    }
#ifdef _DUMP_
    IO::CartesianWriteHDF<float>("result", "S", result, 0);
#endif /* _DUMP_ */

    return 0;
}
