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
