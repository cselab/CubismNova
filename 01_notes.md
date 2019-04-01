<!-- File       : 01_notes.md -->
<!-- Created    : Mon Apr 01 2019 03:06:40 PM (+0200) -->
<!-- Author     : Fabian Wermelinger -->
<!-- Description: Dynamic notebook; remove later -->
<!-- Copyright 2019 ETH Zurich. All Rights Reserved. -->

# Memory management

* Legacy code: sizeX, sizeY and sizeZ is defined at user level -> AoS for
convenience


## Hierarchy (compute node-level)

0. Memory pool [low-level/core]:
  * Memory is viewed as _chunks_ (even if all of data is coalesced)
  * Meta data of a chunk is encoded in a `Node`
  * A `Node` encodes (incomplete) structure with neighbor nodes (connectivity)
  * Allocators live here (POSIX, simple segregated)
  * Depending on the higher level structure, allocation for this structure may
  be carried at once (static objects) or dynamically (tree, MR)
  * Dynamically allocated/extended structures may or may not be coalesced in
  memory (the set of all `Node`s).  A chunk (memory pointed to by a `Node`)
  always is.

1. Structure:
  * This level describes rules that define the structural encoding of `Node`s.
  * Such rules may be dynamic: a local `Node` encoding may change at any later
  time(s).
  * Or static: the local encoding of a `Node` is known prior to memory
  allocation and is guaranteed to remain constant.
  * Defines an allocation scheme/plan
  * This level defines an (dynamic) index space.
  * E.g.: SFC, linked lists, structured (Grid.h/GridMPI.h), trees, multigrid, 
  particles (dynamic structure, how does this map to the `Node` concept (`Node`
  has the same meaning as a cell from a cell-list?)?)

2. Field/Tensor (entry point MPI static allocation plan; dynamic plans here or
   1.[?]):
  * Object with structure, can be indexed into.
  * Defines dimensions (scalar, vector, rank-n tensor in general)
  * Defines data type/precision
  * Defines physical meaning (defines units)
  * Particle fields

3. Block (legacy name) also Fold [high-level/base]:
  * Encodes space
  * Can combine one or multiple tensors
  * At each point (an element of the discrete set space) a tensor represents a
  value
  * Some distinctions:
    - A point may match the tensor value exactly (in index space)
    - If no exact match interpolation is required.
    - Flavors of Blocks: IndexBlock (discrete, more efficient), InterpBlock
