// File       : BlockField.h
// Created    : Thu Apr 11 2019 08:22:16 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Field types with block-structured memory layout
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef BLOCKFIELD_H_6LD9GMJA
#define BLOCKFIELD_H_6LD9GMJA

#include "Core/BlockAllocator.h"
#include "Core/Common.h"

#include <array>
#include <cassert>
#include <cstring>

NAMESPACE_BEGIN()
// String literals used to add a descriptive tag to describe a possible data
// mapping within a block.
//
// Undefined: No particular data layout (default)
// Cell:      Data is represented at cell centers
// Node:      Data is represented at cell nodes (vertices)
// Face:      Data is represented at cell faces (any)
// FaceX:     Data is represented at cell faces in X dimension
// FaceY:     Data is represented at cell faces in Y dimension
// FaceZ:     Data is represented at cell faces in Z dimension
constexpr std::array<const char *, 7> DATA_MAPPING = {
    "Undefined", "Cell", "Node", "Face", "FaceX", "FaceY", "FaceZ"};
NAMESPACE_END()

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(BlockField)

/// @brief Data layout constants used to describe a special block-data layout.
///        Undefined corresponds to no associated mapping (default), Cell
///        corresponds to data mapped to cell centers, Node corresponds to
///        data mapped to cell nodes (vertices), Face corresponds to data
///        mapped to cell faces (independent of direction), FaceX corresponds
///        to data mapped to faces along X dimension, FaceY corresponds to data
///        mapped to faces along Y dimension and FaceZ corresponds to data
///        mapped to faces along Z dimesnion.
enum class DataMapping { Undefined = 0, Cell, Node, Face, FaceX, FaceY, FaceZ };

/// @brief Allocator independent base class for a block with dimensions DIMX,
/// DIMY, DIMZ.
///
/// @tparam DIMX
/// @tparam DIMY
/// @tparam DIMZ
template <size_t DIMX, size_t DIMY, size_t DIMZ>
class BlockBase
{
public:
    static constexpr size_t BlockDimX = DIMX;
    static constexpr size_t BlockDimY = DIMY;
    static constexpr size_t BlockDimZ = DIMZ;

    BlockBase() = default;
    virtual ~BlockBase() {}

    /// @brief Mutable pointer to address of first data element
    ///
    /// @return
    virtual void *getBlockPtr() { return nullptr; }

    /// @brief Immutable pointer to address of first data element
    ///
    /// @return
    virtual const void *getBlockPtr() const { return nullptr; }

    /// @brief Number of bytes occupied by block
    ///
    /// @return Number of bytes
    virtual size_t getBlockBytes() const { return 0; }

    /// @brief Number of data elements contained in block
    ///
    /// @return Number of block data elements
    size_t getBlockSize() const { return DIMX * DIMY * DIMZ; }
};

/// @brief Generic single block field
///
/// @tparam BlockAlloc (block allocator type)
/// @tparam DM         (data mapping identifier)
/// @tparam OWNER      (owner of allocated block memory)
template <typename BlockAlloc,
          DataMapping DM,
          Cubism::Dim DIM,
          bool OWNER = true>
class Field : public BlockBase<BlockAlloc::BlockDimX,
                               BlockAlloc::BlockDimY,
                               BlockAlloc::BlockDimZ>
{
public:
    using BaseType = BlockBase<BlockAlloc::BlockDimX,
                               BlockAlloc::BlockDimY,
                               BlockAlloc::BlockDimZ>;
    using DataType = typename BlockAlloc::DataType;

    static constexpr bool DataOwner = OWNER;
    // Name ID for data element mapping (cell, node, face)
    static constexpr const char *MapName =
        DATA_MAPPING[static_cast<size_t>(DM)];
    // Data mapping class identifier
    static constexpr DataMapping MapClass = DM;
    // Dominant/leading dimension indicator
    static constexpr Cubism::Dim Dim = DIM;

    /// @brief Base constructor
    ///
    /// @param bptr (address to external block memory [only required for
    ///        proxies])
    Field(DataType *bptr = nullptr) : block_(bptr), bytes_(0)
    {
        if (OWNER) {
            // data owner
            assert(
                block_ == nullptr &&
                "Passing a valid address for a field owner is not permitted");
            allocBlock_();
            clearBlock_();
        } else {
            // proxy (nullptr will cause an illegal instruction when
            // dereferenced)
            assert(
                block_ != nullptr &&
                "A field proxy can only be constructed with a valid address");
            bytes_ = blk_alloc_.getBytes(1);
        }
    }

    /// @brief Copy constructor
    ///
    /// @param c
    Field(const Field &c) : block_(nullptr), bytes_(0)
    {
        if (OWNER) {
            // data owner (deep copy)
            allocBlock_();
            copyBlock_(c.block_);
        } else {
            // proxy (shallow copy)
            copyBlockShallow_(c.block_);
        }
    }

    /// @brief Move constructor (move semantics are not permitted for proxies)
    ///
    /// @param c
    Field(Field &&c) noexcept : block_(nullptr), bytes_(0)
    {
        if (OWNER) {
            // move semantics are permitted for data owners
            deallocBlock_();
            copyBlockShallow_(c.block_);
            c.setNull_(); // ensures that destructor has no effect
            return;
        }
        assert(false && "Move semantics in proxy fields are forbidden");
    }

    ~Field()
    {
        if (OWNER) {
            deallocBlock_();
        }
    }

    void *getBlockPtr() override { return static_cast<void *>(block_); }

    const void *getBlockPtr() const override
    {
        return static_cast<const void *>(block_);
    }

    size_t getBlockBytes() const override { return bytes_; }

    /// @brief Copy assignment operator
    ///
    /// @param c
    ///
    /// @return
    Field &operator=(const Field &c)
    {
        if (this != &c) {
            if (OWNER) {
                // data owner (deep copy)
                copyBlock_(c.block_);
            } else {
                // proxy (shallow copy)
                copyBlockShallow_(c.block_);
            }
        }
        return *this;
    }

    /// @brief Move assignment operator (move semantics are not permitted for
    ///        proxies)
    ///
    /// @param c
    ///
    /// @return
    Field &operator=(Field &&c)
    {
        if (OWNER) {
            if (this != &c) {
                // data owner (shallow copy)
                deallocBlock_();
                copyBlockShallow_(c.block_);
                c.setNull_(); // ensures that destructor has no effect
            }
            return *this;
        }
        // a field proxy can never be the owner of data
        assert(false && "Move semantics in proxy fields are forbidden");
        setNull_(); // this will force the proxy to crash when block_ is
                    // dereferenced
        return *this;
    }

    DataType &operator[](size_t i)
    {
        assert(i < this->getBlockSize() && "Linear index out of bounds");
        return block_[i];
    }

    const DataType &operator[](size_t i) const
    {
        assert(i < this->getBlockSize() && "Linear index out of bounds");
        return block_[i];
    }

    DataType &operator()(size_t ix, size_t iy = 0, size_t iz = 0)
    {
        assert(ix < BaseType::BlockDimX && "Block X-index out of bounds");
        assert(iy < BaseType::BlockDimY && "Block Y-index out of bounds");
        assert(iz < BaseType::BlockDimZ && "Block Z-index out of bounds");
        return block_[ix +
                      BaseType::BlockDimX * (iy + BaseType::BlockDimY * iz)];
    }

    const DataType &operator()(size_t ix, size_t iy = 0, size_t iz = 0) const
    {
        assert(ix < BaseType::BlockDimX && "Block X-index out of bounds");
        assert(iy < BaseType::BlockDimY && "Block Y-index out of bounds");
        assert(iz < BaseType::BlockDimZ && "Block Z-index out of bounds");
        return block_[ix +
                      BaseType::BlockDimX * (iy + BaseType::BlockDimY * iz)];
    }

    ////////////////////////////////////////////////////////////////////////////
    // TODO: [fabianw@mavt.ethz.ch; 2019-04-13] arithmetic block operators
    // (outside of class)
    ////////////////////////////////////////////////////////////////////////////

private:
    DataType *block_;      // pointer to first block element
    size_t bytes_;         // number of bytes pointed to by block_
    BlockAlloc blk_alloc_; // block allocator

    /// @brief Allocate single block memory
    ///
    /// @return True if success
    bool allocBlock_()
    {
        block_ = blk_alloc_.allocate(1);
        bytes_ = blk_alloc_.getBytes(1);
        assert(block_ != nullptr && "Block allocator returned NULL address");
        if (block_ == nullptr) {
            return false;
        }
        return true;
    }

    /// @brief Deallocate single block memory
    void deallocBlock_()
    {
        if (block_ != nullptr) {
            blk_alloc_.deallocate(block_);
            bytes_ = 0;
        }
    }

    /// @brief Set all addresses to NULL
    void setNull_()
    {
        block_ = nullptr;
        bytes_ = 0;
    }

    /// @brief Deep copy of an external source
    ///
    /// @param src
    void copyBlock_(const DataType *src)
    {
        if (block_ != src) {
            std::memcpy(block_, src, bytes_);
        }
    }

    /// @brief Shallow copy of an external source (required for move semantics
    ///        and field proxies)
    ///
    /// @param src
    void copyBlockShallow_(DataType *src)
    {
        block_ = src;
        bytes_ = blk_alloc_.getBytes(1);
    }

    /// @brief Clear block memory bitwise
    void clearBlock_()
    {
        if (block_ != nullptr) {
            // DataType should really be built-in or POD
            std::memset(block_, 0, bytes_);
        }
    }
};

/// @brief Cell centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam OWNER (block memory owner)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          bool OWNER = true,
          typename BlockAlloc = AlignedBlockAllocator<DataType, BDX, BDY, BDZ>>
using FieldCell = Field<BlockAlloc, DataMapping::Cell, Cubism::Dim::All, OWNER>;

/// @brief Node centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam OWNER (block memory owner)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          bool OWNER = true,
          typename BlockAlloc =
              AlignedBlockAllocator<DataType, BDX + 1, BDY + 1, BDZ + 1>>
using FieldNode = Field<BlockAlloc, DataMapping::Node, Cubism::Dim::All, OWNER>;

/// @brief Face centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam OWNER (block memory owner)
/// @tparam BlockAlloc (block allocator type)
// template <typename DataType,
//           size_t BDX,
//           size_t BDY,
//           size_t BDZ,
//           bool OWNER = true,
//           typename BlockAlloc = myvec<DataType>>
// // XXX: [fabianw@mavt.ethz.ch; 2019-04-13] Split face fields in X, Y and Z
// using FieldFace = Field<DataType,
//                         DataMapping::Face,
//                         Cubism::Dim::X,
//                         BlockAlloc,
//                         BDX + 1,
//                         BDY + 1,
//                         3 * (BDZ + 1)>;

NAMESPACE_END(BlockField)
NAMESPACE_END(Cubism)

#endif /* BLOCKFIELD_H_6LD9GMJA */
