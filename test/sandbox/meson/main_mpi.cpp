// File       : main_mpi.cpp
// Created    : Tue Mar 23 2021 10:30:43 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Test program using the CubismNova library (MPI)
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include <Cubism/Grid/CartesianMPI.h>
#include <Cubism/Mesh/StructuredUniform.h>
#include <algorithm>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    // 3D mesh
    constexpr size_t dim = 3; // 3D problem
    using Mesh = Cubism::Mesh::StructuredUniform<double, dim>;
    using MIndex = typename Mesh::MultiIndex;

    // cell centered scalar grid using integer data
    constexpr size_t field_rank = 0; // rank-0 field (scalar)
    using Grid = Cubism::Grid::
        CartesianMPI<int, Mesh, Cubism::EntityType::Cell, field_rank>;

    // number of MPI ranks in the topology (2 processes)
    const MIndex nprocs{2, 1, 1};
    const MIndex nblocks(3);     // number of blocks per rank (27 blocks)
    const MIndex block_cells(8); // number of cells per block (512 cells)

    // allocate the grid in the domain [0, 1] (memory is not touched)
    Grid grid(MPI_COMM_WORLD, nprocs, nblocks, block_cells);

    // initialize the block fields using the master thread on each rank with
    // some examples to access state.
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (auto bf : grid) {
        std::fill(bf->begin(), bf->end(), 0); // initialize data to 0
        const auto &fs = bf->getState();      // get state for this field
        const Mesh &bm = *fs.mesh;        // get the block mesh for this field
        const MIndex bi = fs.block_index; // get the block index for this field
        std::cout << "Rank [" << rank << "/" << size
                  << "]: block index: " << bi;
        std::cout << "; cells in block mesh: "
                  << bm.getSize(Cubism::EntityType::Cell) << '\n';
    }

    MPI_Finalize();
    return 0;
}
