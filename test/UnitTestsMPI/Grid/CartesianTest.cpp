// File       : CartesianMPITest.cpp
// Created    : Sun Jan 05 2020 09:27:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian Grid test
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Core/Vector.h"
#include "Grid/CartesianMPI.h"
#include "Mesh/StructuredUniform.h"
#include "gtest-mpi-listener.hpp"
#include "gtest/gtest.h"

#include <mpi.h>

namespace
{
using namespace Cubism;

TEST(CartesianMPI, Construction)
{
    // 3D mesh
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using MIndex = typename Mesh::MultiIndex;
    using Grid = Grid::CartesianMPI<double, Mesh, EntityType::Cell, 0>;
    using PointType = typename Mesh::PointType;

    // number of processes, grid blocks and cells per block
    const MIndex nprocs(2);      // 8 ranks
    const MIndex nblocks(2);     // number of blocks per rank
    const MIndex block_cells(8); // number of cells per block
    Grid grid(MPI_COMM_WORLD, nprocs, nblocks, block_cells);

    const MPI_Comm comm = grid.getCartComm();
    int rank;
    MPI_Comm_rank(comm, &rank);

    EXPECT_EQ(grid.getCartRank(), rank);
    EXPECT_EQ(grid.isRoot(), (0 == rank));
    EXPECT_EQ(grid.getNumProcs(), nprocs);

    const PointType rank_extent = PointType(1) / PointType(nprocs);
    EXPECT_EQ(grid.getMesh().getOrigin(),
              rank_extent * PointType(grid.getProcIndex()));
    EXPECT_EQ(grid.getMesh().getGlobalOrigin(), PointType(0));

    EXPECT_EQ(grid.size(), nblocks.prod());
    EXPECT_EQ(grid.getSize(), nblocks);
    for (auto bf : grid) { // block field in grid
        EXPECT_TRUE(bf->isMemoryOwner());
        EXPECT_NE(bf->getBlockPtr(), nullptr);
        EXPECT_EQ(
            reinterpret_cast<size_t>(bf->getBlockPtr()) % CUBISM_ALIGNMENT, 0);
        EXPECT_NE(&(bf->getState()), nullptr);
        EXPECT_NE((bf->getState()).mesh, nullptr);
    }
}

TEST(CartesianMPI, BlockMesh)
{
    // 2D mesh
    using Mesh = Mesh::StructuredUniform<double, 2>;
    using MIndex = typename Mesh::MultiIndex;
    using Grid = Grid::CartesianMPI<int, Mesh, EntityType::Cell, 0>;

    using PointType = typename Grid::PointType;
    using RealType = typename Grid::RealType;
    using FieldState = typename Grid::FieldState;
    using IntVec = Core::Vector<int, Mesh::Dim>;

    // number of processes, grid blocks and cells per block
    const MIndex nprocs{2, 4};   // 8 ranks
    const MIndex nblocks{2, 5};  // number of blocks per rank
    const MIndex block_cells(8); // number of cells per block
    Grid grid(MPI_COMM_WORLD, nprocs, nblocks, block_cells);

    const MIndex all_blocks = grid.getGlobalSize();
    const MIndex pe_index = grid.getProcIndex();
    EXPECT_EQ(all_blocks, nprocs * nblocks);

    const Mesh &gm = grid.getMesh();
    const PointType O = gm.getOrigin();
    const MIndex Oi = gm.getIndexRange(EntityType::Cell).getBegin();
    const PointType h = gm.getCellSize(0);
    const RealType Vh = gm.getCellVolume(0);
    const PointType block_extent = gm.getExtent() / PointType(nblocks);

    PointType extent(0);
    RealType volume = 0;
    IntVec cells(0);
    IntVec nodes(0);
    std::vector<IntVec> faces(Grid::Dim);
    MIndex blocks(0);

    for (auto bf : grid) { // loop over blocks
        const FieldState &fs = bf->getState();
        const Mesh &fm = *fs.mesh;
        const MIndex global_bi = grid.getGlobalBlockIndex(fs.idx);
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
                const RealType diff =
                    std::fabs((fm.getCellSize(ci) - h).sum() / PointType::Dim);
                EXPECT_LE(diff, std::numeric_limits<RealType>::epsilon());
            }
        }
        if (pe_index[1] == 0 &&
            fs.idx[1] == 0) { // number of global entities along x
            cells[0] += fm.getIndexRange(EntityType::Cell).getExtent()[0];
            nodes[0] += fm.getIndexRange(EntityType::Node).getExtent()[0];
            for (size_t d = 0; d < Grid::Dim; ++d) {
                faces[d][0] +=
                    fm.getIndexRange(EntityType::Face, d).getExtent()[0];
            }
        }
        if (pe_index[0] == 0 &&
            fs.idx[0] == 0) { // number of global entities along y
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
                if (global_bi[d] == all_blocks[d] - 1) {
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

    // global entity counts
    MPI_Comm comm = grid.getCartComm();
    const MIndex gcells = grid.getGlobalSize() * block_cells;
    const MIndex gnodes = gcells + 1;
    std::vector<MIndex> gfaces(Grid::Dim);
    for (size_t d = 0; d < Grid::Dim; ++d) {
        MIndex gface(gcells);
        ++gface[d];
        gfaces[d] = gface;
    }
    for (size_t i = 0; i < Grid::Dim; ++i) {
        int global_count = 0;
        MPI_Allreduce(&cells[i], &global_count, 1, MPI_INT, MPI_SUM, comm);
        EXPECT_EQ(global_count, gcells[i]);

        global_count = 0;
        MPI_Allreduce(&nodes[i], &global_count, 1, MPI_INT, MPI_SUM, comm);
        EXPECT_EQ(global_count, gnodes[i]);

        for (size_t d = 0; d < Grid::Dim; ++d) {
            global_count = 0;
            MPI_Allreduce(
                &faces[d][i], &global_count, 1, MPI_INT, MPI_SUM, comm);
            EXPECT_EQ(global_count, gfaces[d][i]);
        }
    }

    // block counts
    int n = nblocks[0] - 1;
    EXPECT_EQ(blocks[0], nblocks[1] * (n * (n + 1) / 2));
    n = nblocks[1] - 1;
    EXPECT_EQ(blocks[1], nblocks[0] * (n * (n + 1) / 2));
}

} // namespace
