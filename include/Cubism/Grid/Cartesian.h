// File       : Cartesian.h
// Created    : Sun Jan 05 2020 04:34:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian grid composed of block fields
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIAN_H_QBSFTWK7
#define CARTESIAN_H_QBSFTWK7

#include "Alloc/AlignedBlockAllocator.h"
#include "Grid/BlockFieldAssembler.h"
#include <cassert>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)

/// @brief Cartesian block (tensor) field
///
/// @tparam TData Field data type
/// @tparam TMesh Mesh type to be associated with fields
/// @tparam TEntity Cubism::EntityType
/// @tparam RANK Rank of (tensor) fields
/// @tparam TAlloc Allocator for field data
template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity = Cubism::EntityType::Cell,
          size_t RANK = 0,
          template <typename> class TAlloc = AlignedBlockAllocator>
class Cartesian
{
public:
    /// @brief Type of mesh
    using MeshType = TMesh;
    /// @brief Index range type
    using IndexRangeType = typename MeshType::IndexRangeType;
    /// @brief Type for higher dimensional index
    using MultiIndex = typename IndexRangeType::MultiIndex;
    /// @brief Type for physical domain ranges spanned by MeshType
    using RangeType = typename MeshType::RangeType;
    /// @brief Type of point in physical domain (a MeshType::Dim vector of
    ///        float or double)
    using PointType = typename MeshType::PointType;

    /// @brief State (meta data) for individual block fields
    ///
    /// This data structure carries individual meta data information for each
    /// block field in the Cartesian topology.  The mesh pointer points to the
    /// block (sub) mesh if topological information is required.
    struct FieldState {
        size_t rank;
        size_t comp;
        MultiIndex idx;
        MeshType *mesh;
    };

protected:
    /// @brief Type of mesh hull (full mesh or sub-mesh)
    using MeshHull = typename MeshType::MeshHull;
    /// @brief Type of block assembler
    using Assembler =
        BlockFieldAssembler<TEntity, TData, FieldState, MeshType, RANK>;

public:
    /// @brief Block (tensor) field type
    using FieldType = typename Assembler::FieldType;
    /// @brief View type for block fields
    using FieldView = Block::FieldView<FieldType>;
    /// @brief Container type for field views
    using FieldContainer = Block::FieldContainer<FieldView>;

    static constexpr size_t Rank = RANK;
    static constexpr size_t NComponents = FieldType::NComponents;
    static constexpr typename Cubism::EntityType EntityType = TEntity;

    /// @brief Main constructor for a Cartesian block field topology
    ///
    /// @param nblocks Number of blocks
    /// @param block_cells Number of cells in each block
    /// @param start Physical origin for this Cartesian grid (lower left)
    /// @param end Physical end for this Cartesian grid (top right)
    /// @param gorigin Physical global origin
    /// @param block_start Start index of the first block in the topology
    Cartesian(const MultiIndex &nblocks,
              const MultiIndex &block_cells,
              const PointType &start = PointType(0),
              const PointType &end = PointType(1),
              const PointType &gorigin = PointType(0),
              const MultiIndex &block_start = MultiIndex(0))
        : nblocks_(nblocks), block_cells_(block_cells),
          block_range_(block_start, block_start + nblocks)
    {
        // allocate the memory
        alloc_();
        // allocate the global mesh (gorigin and start may be different)
        global_mesh_ = new MeshType(gorigin,
                                    RangeType(start, end),
                                    IndexRangeType(nblocks_ * block_cells_),
                                    MeshHull::FullMesh);
        // assemble the block fields
        assembler_.assemble(data_,
                            *global_mesh_,
                            nblocks_,
                            block_cells_,
                            block_bytes_,
                            component_bytes_);
        // initialize the main field access container
        for (auto tf : assembler_.tensor_fields) {
            fields_.pushBack(new FieldView(*tf));
        }
        // No NUMA touch has been carried out until here.  The user should touch
        // the data based on her/his thread partition strategy in the
        // application.

        assert(fields_.size() == assembler_.field_states.size());
        assert(fields_.size() == assembler_.field_meshes.size());
        assert(fields_.size() == assembler_.tensor_fields.size());
    }

    Cartesian() = delete;
    Cartesian(const Cartesian &c) = delete;
    Cartesian(Cartesian &&c) = default;
    Cartesian &operator=(Cartesian &&c) = delete;

    /// @brief Copy assign field data only (not states or mesh)
    ///
    /// @param c Other Cartesian topology of same type
    Cartesian &operator=(const Cartesian &c)
    {
        assert(size() == c.size());
        for (size_t i = 0; i < c.size(); ++i) {
            fields_[i].copyData(c.fields_[i]);
        }
    }

    virtual ~Cartesian() { dispose_(); }

    // block iterators
    using iterator = typename FieldContainer::iterator;
    using const_iterator = typename FieldContainer::const_iterator;
    using reverse_iterator = typename FieldContainer::reverse_iterator;
    using const_reverse_iterator =
        typename FieldContainer::const_reverse_iterator;

    iterator begin() noexcept { return fields_.begin(); }
    const_iterator begin() const noexcept
    {
        return const_iterator(fields_.begin());
    }
    iterator end() noexcept { return fields_.end(); }
    const_iterator end() const noexcept
    {
        return const_iterator(fields_.end());
    }
    reverse_iterator rbegin() noexcept { return fields_.rbegin(); }
    const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(fields_.rbegin());
    }
    reverse_iterator rend() noexcept { return fields_.rend(); }
    const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(fields_.rend());
    }
    const_iterator cbegin() const noexcept { return fields_.cbegin(); }
    const_iterator cend() const noexcept { return fields_.cend(); }
    const_reverse_iterator crbegin() const noexcept
    {
        return fields_.crbegin();
    }
    const_reverse_iterator crend() const noexcept { return fields_.crend(); }

    /// @brief Returns number of block fields in the topology
    size_t size() const { return fields_.size(); }

    /// @brief Returns size of topology in all dimensions
    MultiIndex getSize() const { return nblocks_; }

    /// @brief Returns the global block index
    ///
    /// @param bi Local Cartesian block index
    MultiIndex getGlobalBlockIndex(const MultiIndex &bi) const
    {
        return block_range_.getBegin() + bi;
    }

    /// @brief Get container of fields
    FieldContainer &getFields() { return fields_; }

    /// @brief Get container of fields
    const FieldContainer &getFields() const { return fields_; }

    /// @brief Get vector of field states
    std::vector<FieldState> &getFieldStates()
    {
        return assembler_.field_states;
    }

    /// @brief Get vector of field states
    const std::vector<FieldState> &getFieldStates() const
    {
        return assembler_.field_states;
    }

    /// @brief Multi-index access fields
    FieldView &operator()(const MultiIndex &p)
    {
        return fields_[block_range_.getFlatIndex(p)];
    }

    /// @brief Multi-index access fields (immutable)
    const FieldView &operator()(const MultiIndex &p) const
    {
        return fields_[block_range_.getFlatIndex(p)];
    }

    /// @brief Linear access fields
    FieldView &operator[](const size_t i)
    {
        assert(i < fields_.size());
        return fields_[i];
    }

    /// @brief Linear access fields (immutable)
    const FieldView &operator[](const size_t i) const
    {
        assert(i < fields_.size());
        return fields_[i];
    }

protected:
    const MultiIndex nblocks_;
    const MultiIndex block_cells_;
    const IndexRangeType block_range_;
    MeshType *global_mesh_;
    FieldContainer fields_;
    Assembler assembler_;

private:
    using BlockData = typename FieldType::BaseType;
    using DataType = typename FieldType::DataType;

    DataType *data_;
    TAlloc<DataType> blk_alloc_;
    size_t block_elements_;
    size_t block_bytes_;
    size_t component_bytes_;
    size_t all_bytes_;

    /// @brief Allocate grid memory
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
        all_bytes_ = nfaces * component_bytes_ * FieldType::NComponents;

        // get the allocation
        assert(all_bytes_ > 0);
        data_ = blk_alloc_.allocate(all_bytes_);
    }

    /// @brief Deallocate grid memory
    void dealloc_() { blk_alloc_.deallocate(data_); }

    void dispose_()
    {
        if (global_mesh_) {
            delete global_mesh_;
        }
        assembler_.dispose();
        for (auto fv : fields_) {
            if (fv) {
                delete fv;
            }
        }
        dealloc_();
    }
};

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity,
          size_t RANK,
          template <typename>
          class TAlloc>
constexpr size_t Cartesian<TData, TMesh, TEntity, RANK, TAlloc>::Rank;

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity,
          size_t RANK,
          template <typename>
          class TAlloc>
constexpr size_t Cartesian<TData, TMesh, TEntity, RANK, TAlloc>::NComponents;

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity,
          size_t RANK,
          template <typename>
          class TAlloc>
constexpr typename Cubism::EntityType
    Cartesian<TData, TMesh, TEntity, RANK, TAlloc>::EntityType;

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* CARTESIAN_H_QBSFTWK7 */
