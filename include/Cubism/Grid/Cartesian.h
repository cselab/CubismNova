// File       : Cartesian.h
// Created    : Sun Jan 05 2020 04:34:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian grid composed of block fields
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIAN_H_QBSFTWK7
#define CARTESIAN_H_QBSFTWK7

#include "Alloc/AlignedBlockAllocator.h"
#include "Common.h"
#include "Grid/BlockFieldAssembler.h"
#include <cassert>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)

/**
 * @brief Cartesian block (tensor) field
 * @tparam T Field data type
 * @tparam Mesh Mesh type to be associated with fields
 * @tparam Entity Entity type
 * @tparam RANK Rank of (tensor) fields
 * @tparam Alloc Allocator for field data
 *
 * @rst
 * Cartesian topology composed of block :ref:`field` for the specified entity
 * type.  As opposed to an individual block :ref:`field`, this class manages a
 * structure of arrays (SoA) memory layout for *all* the blocks in the Cartesian
 * topology instead of just individual blocks.  See the :ref:`cartesianmpi` grid
 * section for a distributed variant of this class.
 * @endrst
 */
template <typename T,
          typename Mesh,
          Cubism::EntityType Entity = Cubism::EntityType::Cell,
          size_t RANK = 0,
          template <typename> class Alloc = AlignedBlockAllocator>
class Cartesian
{
public:
    /** @brief Type of mesh */
    using MeshType = Mesh;
    /** @brief Index range type */
    using IndexRangeType = typename MeshType::IndexRangeType;
    /** @brief Type for higher dimensional index */
    using MultiIndex = typename IndexRangeType::MultiIndex;
    /** @brief Type for physical domain ranges spanned by ``MeshType`` */
    using RangeType = typename MeshType::RangeType;
    /** @brief Type of point in physical domain */
    using PointType = typename MeshType::PointType;
    /** @brief Type float used to describe the mesh topology */
    using RealType = typename MeshType::RealType;

    /**
     * @brief Field state
     *
     * @rst
     * State (meta data) for individual block fields This data structure carries
     * individual meta data information for each block field in the Cartesian
     * topology.  The mesh pointer points to the block (sub) mesh if topological
     * information is required.
     * @endrst
     */
    struct FieldState {
        /** @brief Block index */
        MultiIndex idx;
        /** @brief Block mesh */
        MeshType *mesh;
    };

protected:
    /** @brief Type of mesh hull (full mesh or sub-mesh) */
    using MeshHull = typename MeshType::MeshHull;
    /** @brief Type of block assembler */
    using Assembler =
        BlockFieldAssembler<T, RANK, Entity, FieldState, MeshType>;

public:
    /** @brief Block (scalar, tensor, face) field type */
    using BaseType = typename Assembler::BaseType;
    /** @brief Data type of carried fields */
    using DataType = typename Assembler::DataType;
    /** @brief Container type for field views */
    using FieldContainer = typename Assembler::FieldContainer;

    /** @brief Field dimension */
    static constexpr size_t Dim = MeshType::Dim;
    /** @brief Field rank */
    static constexpr size_t Rank = RANK;
    /** @brief Number of field components */
    static constexpr size_t NComponents = BaseType::NComponents;
    /** @brief Entity type of field */
    static constexpr typename Cubism::EntityType EntityType = Entity;

    /**
     * @brief Default constructor (empty topology)
     */
    Cartesian()
        : nblocks_(0), block_cells_(0), block_range_(0), mesh_(nullptr),
          data_(nullptr)
    {
    }

    /**
     * @brief Main constructor for a Cartesian block field topology
     * @param nblocks Number of blocks
     * @param block_cells Number of cells in each block
     * @param start Physical origin for this Cartesian grid (lower left)
     * @param end Physical end for this Cartesian grid (top right)
     * @param gorigin Physical global origin
     */
    Cartesian(const MultiIndex &nblocks,
              const MultiIndex &block_cells,
              const PointType &start = PointType(0),
              const PointType &end = PointType(1),
              const PointType &gorigin = PointType(0))
        : nblocks_(nblocks), block_cells_(block_cells), block_range_(nblocks),
          mesh_(nullptr), data_(nullptr)
    {
        initTopology_(gorigin, start, end);
    }

    /** @brief Deleted copy constructor */
    Cartesian(const Cartesian &c) = delete;
    /** @brief Deleted move constructor */
    Cartesian(Cartesian &&c) = delete;
    /** @brief Deleted move assignment */
    Cartesian &operator=(Cartesian &&c) = delete;

    /**
     * @brief Copy assign field data only
     * @param c Other Cartesian topology of same type
     *
     * @rst
     * This copies the block data only.
     * @endrst
     */
    Cartesian &operator=(const Cartesian &c)
    {
        assert(size() == c.size());
        for (size_t i = 0; i < c.size(); ++i) {
            assembler_.fields[i].copyData(c.assembler_.fields[i]);
        }
    }

    /** @brief Default destructor */
    virtual ~Cartesian() { dispose_(); }

    /** @brief Block field iterator */
    using iterator = typename FieldContainer::iterator;
    /** @brief Block field iterator */
    using const_iterator = typename FieldContainer::const_iterator;
    /** @brief Reverse block field iterator */
    using reverse_iterator = typename FieldContainer::reverse_iterator;
    /** @brief Reverse block field iterator */
    using const_reverse_iterator =
        typename FieldContainer::const_reverse_iterator;

    /** @return Iterator to first block field */
    iterator begin() noexcept { return assembler_.fields.begin(); }
    /** @return Iterator to first block field */
    const_iterator begin() const noexcept
    {
        return const_iterator(assembler_.fields.begin());
    }
    /** @return Iterator to last block field */
    iterator end() noexcept { return assembler_.fields.end(); }
    /** @return Iterator to last block field */
    const_iterator end() const noexcept
    {
        return const_iterator(assembler_.fields.end());
    }
    /** @return Reverse iterator to first block field */
    reverse_iterator rbegin() noexcept { return assembler_.fields.rbegin(); }
    /** @return Reverse iterator to first block field */
    const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(assembler_.fields.rbegin());
    }
    /** @return Reverse iterator to last block field */
    reverse_iterator rend() noexcept { return assembler_.fields.rend(); }
    /** @return Reverse iterator to last block field */
    const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(assembler_.fields.rend());
    }
    /** @return Iterator to first block field */
    const_iterator cbegin() const noexcept
    {
        return assembler_.fields.cbegin();
    }
    /** @return Iterator to last block field */
    const_iterator cend() const noexcept { return assembler_.fields.cend(); }
    /** @return Reverse iterator to first block field */
    const_reverse_iterator crbegin() const noexcept
    {
        return assembler_.fields.crbegin();
    }
    /** @return Reverse iterator to last block field */
    const_reverse_iterator crend() const noexcept
    {
        return assembler_.fields.crend();
    }

    /**
     * @brief Local size of the grid
     * @return Number of block fields in the local grid
     */
    size_t size() const { return assembler_.fields.size(); }

    /**
     * @brief Local size of the grid in all dimensions
     * @return Number of blocks in all dimensions in the local grid
     */
    MultiIndex getSize() const { return nblocks_; }

    /**
     * @brief Get the global block index
     * @param bi Local Cartesian block index
     * @return Global Cartesian block index
     */
    MultiIndex getGlobalBlockIndex(const MultiIndex &bi) const
    {
        return block_range_.getBegin() + bi;
    }

    /**
     * @return Block index range
     */
    IndexRangeType getBlockRange() const { return block_range_; }

    /**
     * @brief Mesh for the grid
     * @return ``const`` reference to mesh
     *
     * Returns the mesh associated to the Cartesian grid
     */
    const MeshType &getMesh() const
    {
        assert(mesh_ != nullptr);
        return *mesh_;
    }

    /**
     * @brief Field container
     * @return Reference to container of block fields
     *
     * @rst
     * The container has type ``Block::FieldContainer``
     * @endrst
     */
    FieldContainer &getFields() { return assembler_.fields; }

    /**
     * @brief Field container
     * @return ``const`` reference to container of block fields
     *
     * @rst
     * The container has type ``Block::FieldContainer``
     * @endrst
     */
    const FieldContainer &getFields() const { return assembler_.fields; }

    /**
     * @brief Field states
     * @return Reference to vector of field states
     */
    std::vector<FieldState *> &getFieldStates()
    {
        return assembler_.field_states;
    }

    /**
     * @brief Field states
     * @return ``const`` reference to vector of field states
     */
    const std::vector<FieldState *> &getFieldStates() const
    {
        return assembler_.field_states;
    }

    /**
     * @brief Block field access
     * @param p Multi-dimensional block index
     * @return Reference to block field
     */
    BaseType &operator[](const MultiIndex &p)
    {
        assert(assembler_.fields.size() > 0);
        return assembler_.fields[block_range_.getFlatIndex(p)];
    }

    /**
     * @brief Block field access
     * @param p Multi-dimensional block index
     * @return ``const`` reference to block field
     */
    const BaseType &operator[](const MultiIndex &p) const
    {
        assert(assembler_.fields.size() > 0);
        return assembler_.fields[block_range_.getFlatIndex(p)];
    }

    /**
     * @brief Linear block field access
     * @param i One-dimensional block index
     * @return Reference to block field
     */
    BaseType &operator[](const size_t i)
    {
        assert(assembler_.fields.size() > 0);
        assert(i < assembler_.fields.size());
        return assembler_.fields[i];
    }

    /**
     * @brief Linear block field access
     * @param i One-dimensional block index
     * @return ``const`` reference to block field
     */
    const BaseType &operator[](const size_t i) const
    {
        assert(assembler_.fields.size() > 0);
        assert(i < assembler_.fields.size());
        return assembler_.fields[i];
    }

    /**
     * @brief Global size of the grid in all dimensions
     * @return Number of blocks in all dimensions in the global grid
     */
    virtual MultiIndex getGlobalSize() const { return nblocks_; }

protected:
    MultiIndex nblocks_;
    MultiIndex block_cells_;
    IndexRangeType block_range_;

    /**
     * @brief Initialize Cartesian topology
     * @param gorigin Global origin of mesh
     * @param start Lower left point of mesh (rectangular box)
     * @param end Upper right point of mesh (rectangular box)
     * @param nranks Number of ranks in topology
     */
    void initTopology_(const PointType &gorigin,
                       const PointType &start,
                       const PointType &end,
                       const MultiIndex &nranks = MultiIndex(1))
    {
        // allocate the memory
        alloc_();
        // allocate the global mesh (gorigin and start may be different)
        mesh_ = new MeshType(gorigin,
                             RangeType(start, end),
                             IndexRangeType(nblocks_ * block_cells_),
                             MeshHull::FullMesh);
        // assemble the block fields
        assembler_.assemble(data_,
                            *mesh_,
                            block_range_,
                            block_cells_,
                            nranks,
                            block_bytes_,
                            component_bytes_);

        // No NUMA touch has been carried out until here.  The user should touch
        // the data based on her/his thread partition strategy in the
        // application.

        assert(assembler_.fields.size() == assembler_.field_states.size());
        assert(assembler_.fields.size() == assembler_.field_meshes.size());
    }

private:
    using BlockData = typename Assembler::BlockData;

    MeshType *mesh_;
    DataType *data_;
    Assembler assembler_;
    Alloc<DataType> blk_alloc_;
    size_t block_elements_;
    size_t block_bytes_;
    size_t component_bytes_;
    size_t all_bytes_;

    /**
     * @brief Allocate grid memory
     */
    void alloc_()
    {
        if (BlockData::EntityType == Cubism::EntityType::Cell) {
            block_elements_ = block_cells_.prod();
        } else if (BlockData::EntityType == Cubism::EntityType::Node) {
            block_elements_ = (block_cells_ + MultiIndex(1)).prod();
        } else if (BlockData::EntityType == Cubism::EntityType::Face) {
            // XXX: [fabianw@mavt.ethz.ch; 2020-01-05] slightly more than needed
            // but easier for alignment and handling in assembler
            block_elements_ = (block_cells_ + MultiIndex(1)).prod();
        }
        all_bytes_ = block_elements_ * sizeof(DataType);

        // align at CUBISM_ALIGNMENT byte boundary
        all_bytes_ = ((all_bytes_ + CUBISM_ALIGNMENT - 1) / CUBISM_ALIGNMENT) *
                     CUBISM_ALIGNMENT;

        // aligned block bytes (may be larger than the minimum number of bytes
        // needed)
        block_bytes_ = all_bytes_;

        // number of bytes for a single component slice of all blocks in the
        // Cartesian topology
        component_bytes_ = block_bytes_ * nblocks_.prod();

        // if the EntityType of this Cartesian grid is Face, we need to take
        // that into account for the allocated data
        size_t nfaces = 1;
        if (BlockData::EntityType == Cubism::EntityType::Face) {
            nfaces = MeshType::Dim;
        }

        // total number of bytes
        all_bytes_ = nfaces * component_bytes_ * BaseType::NComponents;

        // get the allocation
        assert(all_bytes_ > 0);
        data_ = blk_alloc_.allocate(all_bytes_);
        assert(data_ != nullptr);
    }

    /**
     * @brief Deallocate grid memory
     */
    void dealloc_()
    {
        if (data_) {
            blk_alloc_.deallocate(data_);
        }
    }

    void dispose_()
    {
        assembler_.dispose();
        dealloc_();
        if (mesh_) {
            delete mesh_;
        }
    }
};

template <typename T,
          typename Mesh,
          Cubism::EntityType Entity,
          size_t RANK,
          template <typename>
          class Alloc>
constexpr size_t Cartesian<T, Mesh, Entity, RANK, Alloc>::Dim;

template <typename T,
          typename Mesh,
          Cubism::EntityType Entity,
          size_t RANK,
          template <typename>
          class Alloc>
constexpr size_t Cartesian<T, Mesh, Entity, RANK, Alloc>::Rank;

template <typename T,
          typename Mesh,
          Cubism::EntityType Entity,
          size_t RANK,
          template <typename>
          class Alloc>
constexpr size_t Cartesian<T, Mesh, Entity, RANK, Alloc>::NComponents;

template <typename T,
          typename Mesh,
          Cubism::EntityType Entity,
          size_t RANK,
          template <typename>
          class Alloc>
constexpr typename Cubism::EntityType
    Cartesian<T, Mesh, Entity, RANK, Alloc>::EntityType;

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* CARTESIAN_H_QBSFTWK7 */
