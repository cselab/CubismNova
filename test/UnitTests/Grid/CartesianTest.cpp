// File       : CartesianTest.cpp
// Created    : Sun Jan 05 2020 09:27:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian Grid test
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Grid/Cartesian.h"
#include "Core/Index.h"
#include "Mesh/StructuredUniform.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

TEST(Cartesian, Construction)
{
    // 3D mesh
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using MIndex = typename Mesh::MultiIndex;

    // grid blocks and cells per block
    const MIndex nblocks(2);
    const MIndex block_cells(8);

    { // scalar (rank-0) Cartesian node field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Node, 0>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) { // scalar block field
            EXPECT_TRUE(bf->isMemoryOwner());
            EXPECT_NE(bf->getBlockPtr(), nullptr);
        }
    }

    { // scalar (rank-0) Cartesian face field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Face, 0>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) { // scalar block field
            for (auto d : *bf) { // face direction
                EXPECT_TRUE(d->isMemoryOwner());
                EXPECT_NE(d->getBlockPtr(), nullptr);
            }
        }
    }

    { // vector (rank-1)  Cartesian cell field (double)
        using Grid = Grid::Cartesian<double, Mesh, Cubism::EntityType::Cell, 1>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) {   // tensor block field
            for (auto c : *bf) { // tensor field component
                EXPECT_TRUE(c->isMemoryOwner());
                EXPECT_NE(c->getBlockPtr(), nullptr);
            }
        }
    }

    { // vector (rank-1)  Cartesian face field (float)
        using Grid = Grid::Cartesian<float, Mesh, Cubism::EntityType::Face, 1>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) {   // tensor block field
            for (auto c : *bf) { // tensor field component
                for (auto d : *c) { // face direction
                    EXPECT_TRUE(d->isMemoryOwner());
                    EXPECT_NE(d->getBlockPtr(), nullptr);
                }
            }
        }
    }
}
    const MIndex nblocks(2);
    const MIndex block_cells(8);
    Grid g(nblocks, block_cells);
}

} // namespace
