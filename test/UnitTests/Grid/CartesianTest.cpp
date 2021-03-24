// File       : CartesianTest.cpp
// Created    : Sun Jan 05 2020 09:27:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian Grid test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Block/FieldLab.h"
#include "Cubism/Math.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"
#include <algorithm>
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
            const size_t alignment =
                reinterpret_cast<size_t>(bf->getBlockPtr()) % CUBISM_ALIGNMENT;
            EXPECT_EQ(alignment, 0);
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
                const size_t alignment =
                    reinterpret_cast<size_t>(d->getBlockPtr()) %
                    CUBISM_ALIGNMENT;
                EXPECT_EQ(alignment, 0);
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
                const size_t alignment =
                    reinterpret_cast<size_t>(c->getBlockPtr()) %
                    CUBISM_ALIGNMENT;
                EXPECT_EQ(alignment, 0);
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
        for (auto bf : grid) {      // tensor block field in grid
            for (auto f : *bf) {    // face container component
                EXPECT_NE(&f->getState(), nullptr);
                for (auto c : *f) { // tensor field component of face f
                    EXPECT_TRUE(c->isMemoryOwner());
                    EXPECT_NE(c->getBlockPtr(), nullptr);
                    const size_t alignment =
                        reinterpret_cast<size_t>(c->getBlockPtr()) %
                        CUBISM_ALIGNMENT;
                    EXPECT_EQ(alignment, 0);
                    EXPECT_NE(&(c->getState()), nullptr);
                    EXPECT_NE((c->getState()).mesh, nullptr);
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

TEST(Cartesian, BlockAccess)
{
    // 2D mesh
    using Mesh = Mesh::StructuredUniform<float, 2>;
    using MIndex = typename Mesh::MultiIndex;

    // grid blocks and cells per block
    const MIndex nblocks(2);
    const MIndex block_cells(8);

    { // scalar (rank-0) Cartesian cell field (int)
        using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Cell, 0>;
        Grid grid(nblocks, block_cells);
        auto fields = grid.getIndexFunctor();
        const MIndex f00{0, 0};
        const MIndex f01{0, 1};
        const MIndex f10{1, 0};
        const MIndex f11{1, 1};
        EXPECT_EQ(fields(f00 - f10).getState().block_index, f10);
        EXPECT_EQ(fields(f10 - f10).getState().block_index, f00);
        EXPECT_EQ(fields(f00 + f10).getState().block_index, f10);
        EXPECT_EQ(fields(f10 + f10).getState().block_index, f00);
        EXPECT_EQ(fields(f10 - f01).getState().block_index, f11);
        EXPECT_EQ(fields(f11 - f01).getState().block_index, f10);
        EXPECT_EQ(fields(f11 + f01).getState().block_index, f10);
        EXPECT_EQ(fields(f10 + f01).getState().block_index, f11);
        EXPECT_EQ(fields(f11 + f11).getState().block_index, f00);
        EXPECT_EQ(fields(f00 - f11).getState().block_index, f11);
    }
}

template <Cubism::EntityType Entity, size_t DIM>
void testLab()
{
    static_assert(
        Entity == Cubism::EntityType::Cell ||
            Entity == Cubism::EntityType::Node,
        "Entity must be Cubism::EntityType::Cell or Cubism::EntityType::Node");
    using Mesh = Mesh::StructuredUniform<double, DIM>;
    using Point = typename Mesh::PointType;
    using MIndex = typename Mesh::MultiIndex;
    using Grid = Grid::Cartesian<double, Mesh, Entity, 0>;
    using DataType = typename Grid::DataType;
    using FieldType = typename Grid::BaseType;
    using Lab = Block::FieldLab<FieldType>;
    using Stencil = typename Lab::StencilType;

    // grid blocks and cells per block
    const MIndex nblocks(3);
    const MIndex block_cells(8);
    Grid grid(nblocks, block_cells);
    const Point h = grid.getMesh().getCellSize(0);

    // test function: p \in [0,1]
    auto fexact = [h](const Point &p) -> DataType {
        const DataType fac = 2 * M_PI;
        DataType f = 1;
        int k = 0;
        for (auto x : p) {
            if (Entity == EntityType::Node) {
                if (x < 0) {
                    x += h[k];
                } else if (x > 1) {
                    x -= h[k];
                }
            }
            if (k % 2 == 0) {
                f *= std::sin(fac * x);
            } else {
                f *= std::cos(fac * x);
            }
            ++k;
        }
        return f;
    };

    // initialize grid
    for (auto f : grid) {
        FieldType &bf = *f;
        const Mesh &fm = *bf.getState().mesh;
        for (auto &ci : fm[Entity]) {
            const Point p = fm.getCoords(ci, Entity);
            bf[ci] = fexact(p);
        }
    }

    // setup lab
    Lab lab;
    const Stencil s(-2, 3, true);      // tensorial stencil
    const MIndex sbegin(s.getBegin()); // stencil begin
    const MIndex send(s.getEnd() - 1); // stencil end
    lab.allocate(s, grid[0].getIndexRange());

    // process block fields
    auto findex = grid.getIndexFunctor(); // grid block index functor
    for (auto f : grid) {
        const FieldType &bf = *f;                // reference for field
        lab.loadData(bf.getState().block_index,
                     findex); // load the data lab block

        // create a lab mesh
        const auto lab_range = lab.getActiveLabRange(); // lab index range
        const Mesh &bm = *bf.getState().mesh;           // block mesh
        const MIndex lab_cells =
            bm.getIndexRange(EntityType::Cell).getExtent() + send - sbegin;
        const Mesh mlab(bm.getBegin() + Point(sbegin) * h,
                        bm.getEnd() + Point(send) * h,
                        lab_cells,
                        Mesh::MeshIntegrity::SubMesh);

        // process lab
        for (auto i : lab_range) {
            const Point x = mlab.getCoords(i, Entity);
            const DataType adiff = Cubism::myAbs(lab[i + sbegin] - fexact(x));
            EXPECT_LE(adiff, 8 * std::numeric_limits<DataType>::epsilon());
        }
    }
}

TEST(Cartesian, FieldLab)
{
    testLab<Cubism::EntityType::Cell, 1>();
    testLab<Cubism::EntityType::Cell, 2>();
    testLab<Cubism::EntityType::Cell, 3>();

    testLab<Cubism::EntityType::Node, 1>();
    testLab<Cubism::EntityType::Node, 2>();
    testLab<Cubism::EntityType::Node, 3>();
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
        const PointType O = gm.getBegin();
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
            blocks += fs.block_index;
            EXPECT_TRUE(fm.isSubMesh());
            EXPECT_EQ(fm.getGlobalBegin(), gm.getGlobalBegin());
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
            if (fs.block_index[1] == 0) { // number of global entities along x
                cells[0] += fm.getIndexRange(EntityType::Cell).getExtent()[0];
                nodes[0] += fm.getIndexRange(EntityType::Node).getExtent()[0];
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    faces[d][0] +=
                        fm.getIndexRange(EntityType::Face, d).getExtent()[0];
                }
            }
            if (fs.block_index[0] == 0) { // number of global entities along y
                cells[1] += fm.getIndexRange(EntityType::Cell).getExtent()[1];
                nodes[1] += fm.getIndexRange(EntityType::Node).getExtent()[1];
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    faces[d][1] +=
                        fm.getIndexRange(EntityType::Face, d).getExtent()[1];
                }
            }

            { // block mesh origin
                const PointType mO =
                    O + PointType(fs.block_index) * block_extent;
                const RealType diff =
                    std::fabs((fm.getBegin() - mO).sum() / PointType::Dim);
                EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
            }
            { // block mesh extent
                const RealType diff = std::fabs(
                    (fm.getExtent() - block_extent).sum() / PointType::Dim);
                EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
            }
            { // global index offsets
                const MIndex gc = fm.getIndexRange(EntityType::Cell).getBegin();
                EXPECT_EQ(gc, Oi + fs.block_index * block_cells);
                const MIndex gn = fm.getIndexRange(EntityType::Node).getBegin();
                EXPECT_EQ(gn, Oi + fs.block_index * block_cells);
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    const MIndex gf =
                        fm.getIndexRange(EntityType::Face, d).getBegin();
                    EXPECT_EQ(gf, Oi + fs.block_index * block_cells);
                }
            }
            { // local index extents
                MIndex cell_extent = block_cells;
                MIndex node_extent = cell_extent;
                std::vector<MIndex> face_extents;
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    MIndex face_extent = cell_extent;
                    if (fs.block_index[d] == nblocks[d] - 1) {
                        ++node_extent[d];
                        ++face_extent[d];
                    }
                    face_extents.push_back(face_extent);
                }
                const MIndex ce =
                    fm.getIndexRange(EntityType::Cell).getExtent();
                EXPECT_EQ(ce, cell_extent);
                const MIndex ne =
                    fm.getIndexRange(EntityType::Node).getExtent();
                EXPECT_EQ(ne, node_extent);
                for (size_t d = 0; d < Grid::Dim; ++d) {
                    const MIndex fe =
                        fm.getIndexRange(EntityType::Face, d).getExtent();
                    EXPECT_EQ(fe, face_extents[d]);
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
            const size_t ce = gm.getIndexRange(EntityType::Cell).getExtent()[i];
            EXPECT_EQ(ce, cells[i]);
            const size_t ne = gm.getIndexRange(EntityType::Node).getExtent()[i];
            EXPECT_EQ(ne, nodes[i]);
            for (size_t d = 0; d < Grid::Dim; ++d) {
                const size_t fe =
                    gm.getIndexRange(EntityType::Face, d).getExtent()[i];
                EXPECT_EQ(fe, faces[d][i]);
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
