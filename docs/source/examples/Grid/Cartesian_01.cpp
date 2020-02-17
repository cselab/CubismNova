#include "Grid/Cartesian.h"
#include "Mesh/StructuredUniform.h"
#include <algorithm>

using namespace Cubism;

int main(void)
{
    // 3D mesh
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using MIndex = typename Mesh::MultiIndex;

    // cell centered scalar grid using integer data
    using Grid = Grid::Cartesian<int, Mesh, EntityType::Cell, 0>;

    const MIndex nblocks(3);     // number of blocks in the topology (27 blocks)
    const MIndex block_cells(8); // number of cells per block (512 cells)

    // allocate the grid in the domain [0, 1] (memory is not touched)
    Grid grid(nblocks, block_cells);

    // initialize the block fields using the master thread with some examples to
    // access state.
    for (auto bf : grid) {
        std::fill(bf->begin(), bf->end(), 0); // initialize data to 0
        const auto &fs = bf->getState();      // get state for this field
        const Mesh &bm = *fs.mesh; // get the block mesh for this field
        const MIndex bi = fs.block_index; // get the block index for this field
    }

    return 0;
}
