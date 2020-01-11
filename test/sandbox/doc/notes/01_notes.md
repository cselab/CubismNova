<!-- File       : 01_notes.md -->
<!-- Created    : Mon Apr 01 2019 03:06:40 PM (+0200) -->
<!-- Author     : Fabian Wermelinger -->
<!-- Description: Dynamic notebook; remove later -->
<!-- Copyright 2019 ETH Zurich. All Rights Reserved. -->

<!-- vim-markdown-toc GFM -->

* [Thoughts](#thoughts)
  * [Mesh](#mesh)
    * [Interface components](#interface-components)
* [Memory management (low-level thoughts)](#memory-management-low-level-thoughts)
  * [Hierarchy (compute node-level)](#hierarchy-compute-node-level)
  * [Plans](#plans)
  * [Field/Tensor](#fieldtensor)

<!-- vim-markdown-toc -->

# Thoughts

* How are Face/Node fields handled in AMR? (are boundary elements treated as
duplicates?)
* Single indexer/range iterators for any field
* Field data structures must support vectorization
* Block compounds -> variant or dynamic cast (not favorable)
* Geometry involves most of index conversion cell -> face -> node etc.
* Geometry should it inherit from field or wrap
* Geometry data compute on the fly or store in memory?

## Mesh

### Interface components

* getCoords: For entity, flat index, multi index or iterator
* getCellVolume: cell index or iterator
* getCellWidth: cell index or iterator
* getFaceArea: face index, Dir:: or iterator
* getFaceNormal: face index, Dir:: or iterator
* getNeighbors: For entity from and index to entity or from iterator to entity

# Memory management (low-level thoughts)

* Legacy code: sizeX, sizeY and sizeZ is defined at user level -> AoS for
convenience
* Here: fields define the underlying (initial) block structure, type and size
* Memory nodes (pages) are the elementary allocation element (they define no
type, just bytes)
* Pages are castable to type which gives them meaning and dimension
* Pages may be suitable for data compression


## Hierarchy (compute node-level)

0. Memory pool [low-level/core]:
  * Memory is viewed as _chunks_/_page_/_node_ (even if all of data is coalesced)
  * Meta data for a chunk is encoded in a POD type `Page`
  * A `Page` encodes (incomplete) structure with neighbor nodes (connectivity)
  in more complex structures.
  * Allocators live here (POSIX, simple segregated)
  * Depending on the higher level structure, allocation for this structure may
  be carried out at once (static objects) or dynamically (tree, MR)
  * Dynamically allocated/extended structures may or may not be coalesced in
  memory (the set of all `Page`s).  A page always is.

1. Structure:
  * This level describes rules that define the structural encoding of `Page`s.
  * Such rules may be dynamic: a local `Page` encoding may change at any later
  time(s).
  * Or static: the local encoding of a `Page` is known prior to memory
  allocation and is guaranteed to remain constant.
  * Defines an allocation scheme/plan
  * This level defines an (dynamic) index space.
  * E.g.: SFC, linked lists, structured (Grid.h/GridMPI.h), trees, multigrid, 
  particles (dynamic structure, how does this map to the `Page` concept (`Page`
  has the same meaning as a cell from a cell-list?)?)

2. Field/Tensor (entry point MPI static allocation plan; dynamic plans here or
   1.[?]):
  * Object with structure, can be indexed into pages and assign meaning to them.
  * Defines size (dimensions: scalar, vector, rank-n tensor in general)
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
  * special blocks may be emitted from different fields: E.g.
    Only `FieldCartesianStretched` types can emit stretched blocks (space
    encoding is different).  This would possibly allow to have uniform
    (non-stretched) and stretched fields at the same time.


## Plans

`Plan`s are low-level representations of various allocation schemes.  These
schemes may be static or dynamic and operate on a smallest element called
`Page`.  The page defines a block-structured memory layout and is a _linear_
entity, therefore has structure of an array.  The static/dynamic property is
the foundation for AMR.  Plans with simple segregated allocators may be an
option for such a case concerning efficiency.  A plan may manage multiple pages
side-by-side to add dimensionality.  Its basic type may look like:

```
template <size_t Dim, size_t PageByte>
class PlanBasic

template <size_t Dim, size_t PageByte>
class PlanMultigrid

template <size_t Dim, size_t PageByte>
class PlanAMR

template <size_t Dim, size_t PageByte>
class PlanGPU

template <size_t Dim, size_t PageByte>
class PlanHetreogeneous

template <size_t Dim, size_t PageByte>
class PlanParticles
```

etc.  I don't see much fit for a common interface here.  POD for `Page` may
look different for each of these plans.


## Field/Tensor

A field uses a certain `Plan` to accomplish its target.  A common interface may
be defined by

```
template <size_t Dim, typename T>
Field
```

Common interface may include:
* Type access
* Units
* Dimension (scalar [Dim=1], vector [Dim=3], rank-n in general [Dim=n])
* Some iterator for primitive block type (could be used as a generic way for
 dump operations)

A very first derived field would be Cartesian (closest to legacy code)
```
template <size_t Dim, typename T, some enum for StorageType, size_t BSX, size_t BSY=BSX, size_t BSZ=BSX>
FieldCartesian : public Field<Dim, T>
```
where `StorageType` specifies the data representation (cell, face, vertex).
Such a field may use different plans  if complexity/concept
allows.  Perfect polymorphism is probably not possible/difficult/too complex to
cover all possibilities, different field types will be necessary most likely.
