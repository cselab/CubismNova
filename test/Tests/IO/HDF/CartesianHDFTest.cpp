// File       : CartesianHDFTest.cpp
// Created    : Wed Jan 29 2020 10:53:13 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian grid HDF IO
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/IO/CartesianHDF.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"
#include <algorithm>

namespace
{
using namespace Cubism;

TEST(IO, CartesianWriteHDF)
{
    {
        // 2D mesh
        using Mesh = Mesh::StructuredUniform<double, 2>;
        using MIndex = typename Mesh::MultiIndex;
        using PointType = typename Mesh::PointType;

        // grid blocks and cells per block
        const MIndex nblocks(3);
        const MIndex block_cells(8);

        // Scalar (rank-0) Cartesian cell field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Cell, 0>;
        using DataType = typename Grid::DataType;
        Grid grid(nblocks, block_cells);

        DataType k = 0;
        for (auto bf : grid) { // fill scalar block fields with data
            std::fill(bf->begin(), bf->end(), k++);
        }

        // all of the grid
        IO::CartesianWriteHDF<DataType>("g2D", "block_id", grid);

        // sub-mesh selections
        const auto &m = grid.getMesh();
        const auto block = m.getSubMesh(PointType(0.3), PointType(0.7));
        const auto slice = m.getSubMesh(PointType(0.3), PointType{0.3, 0.7});
        const auto xxl = m.getSubMesh(PointType(-10.0), PointType(10.0));
        const auto lefty =
            m.getSubMesh(PointType{-10.0, 0.0}, PointType{10.0, 0.0});
        const auto righty =
            m.getSubMesh(PointType{-10.0, 1.0}, PointType{10.0, 1.0});
        const auto linez =
            m.getSubMesh(PointType{0.5, -10.0}, PointType{0.5, 10.0});
        IO::CartesianWriteHDF<DataType>("g2D_block", "apples", grid, *block);
        IO::CartesianWriteHDF<DataType>("g2D_slice", "bananas", grid, *slice);
        IO::CartesianWriteHDF<DataType>("g2D_xxl", "peaches", grid, *xxl);
        IO::CartesianWriteHDF<DataType>("g2D_lefty", "oranges", grid, *lefty);
        IO::CartesianWriteHDF<DataType>("g2D_righty", "kiwis", grid, *righty);
        IO::CartesianWriteHDF<DataType>("g2D_linez", "beers", grid, *linez);
    }
    {
        // 3D mesh
        using Mesh = Mesh::StructuredUniform<double, 3>;
        using MIndex = typename Mesh::MultiIndex;
        using PointType = typename Mesh::PointType;

        // grid blocks and cells per block
        const MIndex nblocks(3);
        const MIndex block_cells(8);

        // Vector (rank-1) Cartesian node field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Node, 1>;
        using DataType = typename Grid::DataType;
        Grid grid(nblocks, block_cells);

        DataType k = 0;
        for (auto bf : grid) { // fill scalar block fields with data
            for (auto c : *bf) {
                std::fill(c->begin(), c->end(), k++);
            }
        }

        // all of the grid
        IO::CartesianWriteHDF<DataType>("gv3D", "block_id", grid);

        // sub-mesh selections
        const auto &m = grid.getMesh();
        const auto block = m.getSubMesh(PointType(0.3), PointType(0.7));
        const auto slice =
            m.getSubMesh(PointType(0.3), PointType{0.3, 0.7, 0.7});
        const auto xxl = m.getSubMesh(PointType(-10.0), PointType(10.0));
        const auto lefty = m.getSubMesh(PointType{-10.0, 0.0, -10.0},
                                        PointType{10.0, 0.0, 10.0});
        const auto righty = m.getSubMesh(PointType{-10.0, 1.0, -10.0},
                                         PointType{10.0, 1.0, 10.0});
        const auto linez =
            m.getSubMesh(PointType{0.5, 0.5, -10.0}, PointType{0.5, 0.5, 10.0});
        IO::CartesianWriteHDF<DataType>("gv3D_block", "apples", grid, *block);
        IO::CartesianWriteHDF<DataType>("gv3D_slice", "bananas", grid, *slice);
        IO::CartesianWriteHDF<DataType>("gv3D_xxl", "peaches", grid, *xxl);
        IO::CartesianWriteHDF<DataType>("gv3D_lefty", "oranges", grid, *lefty);
        IO::CartesianWriteHDF<DataType>("gv3D_righty", "kiwis", grid, *righty);
        IO::CartesianWriteHDF<DataType>("gv3D_linez", "beers", grid, *linez);
    }
}
} // namespace
