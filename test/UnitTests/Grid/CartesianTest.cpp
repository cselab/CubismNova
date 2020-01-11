// File       : CartesianTest.cpp
// Created    : Sun Jan 05 2020 09:27:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian Grid test
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Grid/Cartesian.h"
#include "Mesh/StructuredUniform.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

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
            EXPECT_NE(&(bf->getState()), nullptr);
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
                EXPECT_NE(&(d->getState()), nullptr);
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
                EXPECT_NE(&(c->getState()), nullptr);
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
                    EXPECT_NE(&(d->getState()), nullptr);
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

        // check
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

TEST(Cartesian, BlockMesh)
{
    // 2D mesh
    using Mesh = Mesh::StructuredUniform<float, 2>;
    using MIndex = typename Mesh::MultiIndex;

    // grid blocks and cells per block
    const MIndex nblocks{4, 7};
    const MIndex block_cells(8);

    { // scalar (rank-0) Cartesian cell field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Cell, 0>;
        using PointType = typename Grid::PointType;
        using RealType = typename Grid::RealType;
        using FieldState = typename Grid::FieldState;
        Grid grid(nblocks, block_cells);
        const Mesh &gm = grid.getMesh();
        const PointType O = gm.getOrigin();
        const MIndex Oi = gm.getIndexRange(EntityType::Cell).getBegin();
        const PointType h = gm.getCellSize(0);
        const RealType Vh = gm.getCellVolume(0);
        const PointType block_extent = gm.getExtent() / PointType(nblocks);
        PointType extent(0);
        RealType volume = 0;
        MIndex cells(0);
        MIndex nodes(0);
        std::vector<MIndex> faces(Grid::Dim);
        MIndex blocks(0);
        for (auto bf : grid) { // loop over blocks
            const FieldState &fs = bf->getState();
            const Mesh &fm = *fs.mesh;
            extent += fm.getExtent();
            volume += fm.getVolume();
            blocks += fs.idx;
            EXPECT_TRUE(fm.isSubMesh());
            EXPECT_EQ(fm.getGlobalOrigin(), gm.getGlobalOrigin());
            for (const auto &ci : fm[EntityType::Cell]) { // cell checks
                {
                    const RealType diff = std::fabs(fm.getCellVolume(ci) - Vh);
                    EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
                }
                {
                    const RealType diff = std::fabs(
                        (fm.getCellSize(ci) - h).sum() / PointType::Dim);
                    EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
                }
            }
            if (fs.idx[1] == 0) { // number of global entities along x
                cells[0] += fm.getIndexRange(EntityType::Cell).getExtent()[0];
                nodes[0] += fm.getIndexRange(EntityType::Node).getExtent()[0];
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    faces[d][0] +=
                        fm.getIndexRange(EntityType::Face, d).getExtent()[0];
                }
            }
            if (fs.idx[0] == 0) { // number of global entities along y
                cells[1] += fm.getIndexRange(EntityType::Cell).getExtent()[1];
                nodes[1] += fm.getIndexRange(EntityType::Node).getExtent()[1];
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    faces[d][1] +=
                        fm.getIndexRange(EntityType::Face, d).getExtent()[1];
                }
            }

            { // block mesh origin
                const PointType mO = O + PointType(fs.idx) * block_extent;
                const RealType diff =
                    std::fabs((fm.getOrigin() - mO).sum() / PointType::Dim);
                EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
            }
            { // block mesh extent
                const RealType diff = std::fabs(
                    (fm.getExtent() - block_extent).sum() / PointType::Dim);
                EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
            }
            { // global index offsets
                EXPECT_EQ(fm.getIndexRange(EntityType::Cell).getBegin(),
                          Oi + fs.idx * block_cells);
                EXPECT_EQ(fm.getIndexRange(EntityType::Node).getBegin(),
                          Oi + fs.idx * block_cells);
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    EXPECT_EQ(fm.getIndexRange(EntityType::Face, d).getBegin(),
                              Oi + fs.idx * block_cells);
                }
            }
            { // local index extents
                MIndex cell_extent = block_cells;
                MIndex node_extent = cell_extent;
                std::vector<MIndex> face_extents;
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    MIndex face_extent = cell_extent;
                    if (fs.idx[d] == nblocks[d] - 1) {
                        ++node_extent[d];
                        ++face_extent[d];
                    }
                    face_extents.push_back(face_extent);
                }
                EXPECT_EQ(fm.getIndexRange(EntityType::Cell).getExtent(),
                          cell_extent);
                EXPECT_EQ(fm.getIndexRange(EntityType::Node).getExtent(),
                          node_extent);
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    EXPECT_EQ(fm.getIndexRange(EntityType::Face, d).getExtent(),
                              face_extents[d]);
                }
            }
        }
        extent /= PointType{nblocks[1], nblocks[0]};
        { // global extent
            const RealType diff =
                std::fabs((extent - gm.getExtent()).sum() / PointType::Dim);
            EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
        }
        { // global volume
            const RealType diff = std::fabs(volume - gm.getVolume());
            EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
        }

        // entity counts
        for (size_t i = 0; i < Grid::Dim; ++i) {
            EXPECT_EQ(cells[i],
                      gm.getIndexRange(EntityType::Cell).getExtent()[i]);
            EXPECT_EQ(nodes[i],
                      gm.getIndexRange(EntityType::Node).getExtent()[i]);
            for (size_t d = 0; d < Grid::Dim; ++d) {
                EXPECT_EQ(faces[d][i],
                          gm.getIndexRange(EntityType::Face, d).getExtent()[i]);
            }
        }

        // block counts
        int n = nblocks[0] - 1;
        EXPECT_EQ(blocks[0], nblocks[1] * (n * (n + 1) / 2));
        n = nblocks[1] - 1;
        EXPECT_EQ(blocks[1], nblocks[0] * (n * (n + 1) / 2));
    }
}

} // namespace
