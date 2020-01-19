#include "Grid/CartesianMPI.h"
#include "Mesh/StructuredUniform.h"
#include <algorithm>
#include <mpi.h>

using namespace Cubism;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    // 3D mesh
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using MIndex = typename Mesh::MultiIndex;

    // cell centered scalar grid using integer data
    using Grid = Grid::CartesianMPI<int, Mesh, EntityType::Cell, 0>;

    const MIndex nprocs(2); // number of MPI ranks in the topology (8 processes)
    const MIndex nblocks(3);     // number of blocks per rank (27 blocks)
    const MIndex block_cells(8); // number of cells per block (512 cells)

    // allocate the grid in the domain [0, 1] (memory is not touched)
    Grid grid(MPI_COMM_WORLD, nprocs, nblocks, block_cells);

    // initialize the block fields using the master thread on each rank
    for (auto bf : grid) {
        std::fill(bf->begin(), bf->end(), 0);
    }

    MPI_Finalize();
    return 0;
}
