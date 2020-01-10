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
    const MIndex nblocks(3);
    const MIndex block_cells(8);

    { // scalar (rank-0) Cartesian node field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Node, 0>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) { // scalar block field in grid
            EXPECT_TRUE(bf->isMemoryOwner());
            EXPECT_NE(bf->getBlockPtr(), nullptr);
            EXPECT_EQ(reinterpret_cast<size_t>(bf->getBlockPtr()) %
                          CUBISM_ALIGNMENT,
                      0);
            EXPECT_NE((bf->getState()).mesh, nullptr);
        }
    }

    { // scalar (rank-0) Cartesian face field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Face, 0>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) {   // scalar block field in grid
            for (auto d : *bf) { // face direction
                EXPECT_TRUE(d->isMemoryOwner());
                EXPECT_NE(d->getBlockPtr(), nullptr);
                EXPECT_EQ(reinterpret_cast<size_t>(d->getBlockPtr()) %
                              CUBISM_ALIGNMENT,
                          0);
                EXPECT_NE((d->getState()).mesh, nullptr);
            }
        }
    }

    { // vector (rank-1)  Cartesian cell field (double)
        using Grid = Grid::Cartesian<double, Mesh, Cubism::EntityType::Cell, 1>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) {   // tensor block field in grid
            for (auto c : *bf) { // tensor field component
                EXPECT_TRUE(c->isMemoryOwner());
                EXPECT_NE(c->getBlockPtr(), nullptr);
                EXPECT_EQ(reinterpret_cast<size_t>(c->getBlockPtr()) %
                              CUBISM_ALIGNMENT,
                          0);
                EXPECT_NE((c->getState()).mesh, nullptr);
            }
        }
    }

    { // vector (rank-1)  Cartesian face field (float)
        using Grid = Grid::Cartesian<float, Mesh, Cubism::EntityType::Face, 1>;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        for (auto bf : grid) {   // tensor block field in grid
            for (auto c : *bf) { // tensor field component
                for (auto d : *c) { // face direction
                    EXPECT_TRUE(d->isMemoryOwner());
                    EXPECT_NE(d->getBlockPtr(), nullptr);
                    EXPECT_EQ(reinterpret_cast<size_t>(d->getBlockPtr()) %
                                  CUBISM_ALIGNMENT,
                              0);
                    EXPECT_NE((d->getState()).mesh, nullptr);
                }
            }
        }
    }
}

TEST(Cartesian, GridFill)
{
    // 2D mesh
    using Mesh = Mesh::StructuredUniform<float, 2>;
    using MIndex = typename Mesh::MultiIndex;

    // grid blocks and cells per block
    const MIndex nblocks(2);
    const MIndex block_cells(8);

    { // scalar (rank-0) Cartesian cell field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Cell, 0>;
        using DataType = typename Grid::DataType;
        Grid grid(nblocks, block_cells);
        EXPECT_EQ(grid.size(), nblocks.prod());
        EXPECT_EQ(grid.getSize(), nblocks);
        DataType k = 0;
        for (auto bf : grid) { // fill scalar block fields with data
            std::fill(bf->begin(), bf->end(), k++);
        }

        k = nblocks.prod() - 1;
        const DataType ref = block_cells.prod() * k * (k + 1) / 2;
        k = 0;
        for (auto bf : grid) {   // sum all cell values
            for (auto c : *bf) { // for each value in block field bf
                k += c;
            }
        }
        EXPECT_EQ(k, ref);
    }
}
}

} // namespace
