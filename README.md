Refactored code based on the [Cubism](https://gitlab.ethz.ch/mavt-cse/Cubism.git)
ancestor.

# Major changes

* Flexible representation of data with support for space filling curves
multigrid and tree structures for AMR.

* Introduction of fields for the representation of multiple data that map to
the same spatial coordinates.

* Fields allow for cell centered, face centered or cell vertex representations
of data.

* Fields are stored in memory as structures of arrays (SoA) for better
separation and organization in container objects (blocks [rectangular cuboid])
that possibly combine multiple fields for operations on the data.

* Enhanced interface for mesh geometry and indexing that maps to the
corresponding data representation of the field.

* Python bindings for easier configuration and user interaction efficiency.

* Documentation
  * [readthedocs](https://readthedocs.org/)
  * [doxygen](http://www.doxygen.org/)


# Additional differences

* Support for Co-routines (yield, affects MPI communication mainly, recursion)

* Provide a set of basic kernels to perform operations on data:
  * Finite difference schemes
  * Prolongation/restriction for multigrid
  * Interpolation between cell centered, face centered and cell vertex
  representations
  * Volume rendering (?)
  * Marching cubes (?)

* Support for ISPC, GPU (CUDA) and MATE (code generation in general[?])

# Drawbacks of the old version

* Rigid communication sequence mirrored both in the kernels and the 
block processor.

  Consider an algorithm that requires communication
  (e.g. exchange of halo cells, reduction, linear solver, dump).
  Each MPI rank operates on multiple blocks
  and executes the same kernel code on many blocks.
  A kernel cannot issue calls to MPI or other libraries because

  - the call would be performed multiple times;
  - a block is not aware of other blocks so it cannot coordinate communication.

  Therefore, a communication request needs to be reflected both
  in the kernels (they require communication as prerequisites)
  and the block processor (they ensure communication before executing the kernels).

  Drawbacks

  - maintaining two separate parts of the implementation:
  kernels and block processor
  - limited code reuse
  (e.g. algorithms that require communication cannot be
  encapsulated in a function)

* All kernels share a common data structure.

  Particularly applies to MPCF with its FluidElement. 
  The kernels (e.g. HLLC and viscous fluxes) together with their counterparts
  on the node level potentially have access to all fields in the FluidElement.
  Any temporary field also needs to be included in the FluidElement.
  
  Drawbacks

  - kernels assume the presence of certain fields and
  cannot be reused in a different context 
  - temporary fields required by the kernels are also shared 
  which complicates their interaction
  - no well-defined interfaces to kernels as all data is exposed

* No support for indexing and mesh topology.

  Kernels on the core level effectively operate on a 3D array
  taking as arguments the dimensions and a pointer to data.
  
  Drawbacks

  - stencil operations are only expressed using indices such as `[i][j][k]` 
  which is cumbersome and error prone.
