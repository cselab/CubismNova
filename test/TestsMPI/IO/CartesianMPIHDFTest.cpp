// File       : CartesianMPIHDFTest.cpp
// Created    : Tue Feb 04 2020 12:44:45 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian MPI grid HDF IO
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/IO/CartesianMPIHDF.h"
#include "Cubism/Grid/CartesianMPI.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <mpi.h>

namespace
{
using namespace Cubism;

template <Cubism::EntityType Entity>
struct Initializer {
    template <typename Grid>
    void init(Grid &g, const typename Grid::DataType o)
    {
        typename Grid::DataType k = 0;
        for (auto f : g) {
            std::fill(f->begin(), f->end(), o + k++);
        }
    }
};

template <>
struct Initializer<Cubism::EntityType::Face> {
    template <typename Grid>
    void init(Grid &g, const typename Grid::DataType o)
    {
        typename Grid::DataType k = 0;
        for (auto ff : g) {      // face container
            for (auto f : *ff) { // face direction
                std::fill(f->begin(), f->end(), o + k++);
            }
        }
    }
};

TEST(IO, CartesianMPIWriteHDF)
{
    {
        // 2D mesh
        using Mesh = Mesh::StructuredUniform<double, 2>;
        using MIndex = typename Mesh::MultiIndex;
        using PointType = typename Mesh::PointType;

        using DataType = int;
        using Grid =
            Grid::CartesianMPI<DataType, Mesh, Cubism::EntityType::Cell, 0>;

        // grid blocks and cells per block
        const MIndex nprocs{4, 2}; // 8 ranks
        const MIndex nblocks(3);
        const MIndex block_cells(8);

        Grid grid(MPI_COMM_WORLD, nprocs, nblocks, block_cells);
        const MPI_Comm comm = grid.getCartComm();
        int rank;
        MPI_Comm_rank(comm, &rank);
        Initializer<Grid::EntityType> ginit;
        ginit.init(grid, nblocks.prod() * rank);

        // all of the grid
        const auto &m = grid.getGlobalMesh();

        // sub-mesh selections
        //
        // sub-block (from 0.3 to 0.7)
        const auto blk = m.getSubMesh(PointType(0.3), PointType(0.7));
        // slice at y=0 (left boundary)
        const auto lfy =
            m.getSubMesh(PointType{-10.0, 0.0}, PointType{10.0, 0.0});

        const double time = 0;
        IO::CartesianMPIWriteHDF<DataType>("cg2Dmpiall", "basket", grid, time);
        IO::CartesianMPIWriteHDF<DataType>(
            "cg2Dmpiblk", "apples", grid, *blk, time);
        IO::CartesianMPIWriteHDF<DataType>(
            "cg2Dmpilfy", "oranges", grid, *lfy, time);
    }
}

TEST(IO, CartesianWriteReadBackHDF)
{
    // 3D mesh
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using MIndex = typename Mesh::MultiIndex;
    using DataType = int;
    using Grid =
        Grid::CartesianMPI<DataType, Mesh, Cubism::EntityType::Cell, 0>;

    // grid blocks and cells per block
    const MIndex nprocs(2); // 8 ranks
    const MIndex nblocks(3);
    const MIndex block_cells(8);

    Grid src(MPI_COMM_WORLD, nprocs, nblocks, block_cells);
    const MPI_Comm comm = src.getCartComm();
    int rank;
    MPI_Comm_rank(comm, &rank);
    Initializer<Grid::EntityType> ginit;
    const int my_offset = nblocks.prod() * rank;
    ginit.init(src, my_offset);

    // write the field
    IO::CartesianMPIWriteHDF<DataType>("framenmpi", "ramen", src, 0);

    // read back soup
    Grid dst(MPI_COMM_WORLD, nprocs, nblocks, block_cells);
    IO::CartesianMPIReadHDF<DataType>("framenmpi", dst);

    int k = 0;
    for (const auto cf : dst) {
        for (const auto c : *cf) {
            EXPECT_EQ(c, my_offset + k);
        }
        ++k;
    }
}
} // namespace
