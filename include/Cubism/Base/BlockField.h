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
#include <vector>

NAMESPACE_BEGIN()
// String literals used to add a descriptive tag to describe a possible data
// mapping within a block.
//
// Undefined: No particular data layout (default)
// Cell:      Data is represented at cell centers
// Node:      Data is represented at cell nodes (vertices)
// Face:      Data is represented at cell faces (any)
constexpr std::array<const char *, 4> DATA_MAPPING = {
    "Undefined", "Cell", "Node", "Face"};
NAMESPACE_END()

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(BlockField)

/// @brief Data layout constants used to describe a special block-data layout.
///        Undefined corresponds to no associated mapping (default), Cell
///        corresponds to data mapped to cell centers, Node corresponds to
///        data mapped to cell nodes (vertices), Face corresponds to data
///        mapped to cell faces
enum class DataMapping { Undefined = 0, Cell, Node, Face };

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

/// @brief Generic single block field that owns the block memory
///
/// @tparam BlockAlloc (block allocator type)
/// @tparam DM         (data mapping identifier)
/// @tparam DIR        (direction indicator)
template <typename BlockAlloc, DataMapping DM, Cubism::Dir DIR>
class Field : public BlockBase<BlockAlloc::BlockDimX,
                               BlockAlloc::BlockDimY,
                               BlockAlloc::BlockDimZ>
{
public:
    using BaseType = BlockBase<BlockAlloc::BlockDimX,
                               BlockAlloc::BlockDimY,
                               BlockAlloc::BlockDimZ>;
    using AllocType = BlockAlloc;
    using DataType = typename BlockAlloc::DataType;

    // Name ID for data element mapping (cell, node, face)
    static constexpr const char *MapName =
        DATA_MAPPING[static_cast<size_t>(DM)];
    static constexpr DataMapping MapClass = DM; // Data mapping class identifier
    static constexpr Cubism::Dir Dir = DIR;     // Direction indicator

    /// @brief Base constructor
    Field(const bool alloc = true) : block_(nullptr), bytes_(0)
    {
        if (alloc) {
            assert(
                block_ == nullptr &&
                "Passing a valid address for a field owner is not permitted");
            allocBlock_();
            // TODO: [fabianw@mavt.ethz.ch; 2019-04-14] Works only with POD
            // currently
            clearBlock_();
        }
    }

    /// @brief Copy constructor
    Field(const Field &c) : block_(nullptr), bytes_(0)
    {
        allocBlock_();
        copyBlock_(c.block_); // deep copy
    }

    /// @brief Move constructor
    Field(Field &&c) noexcept : block_(nullptr), bytes_(0)
    {
        deallocBlock_();
        copyBlockShallow_(c.block_);
        c.setNull_(); // ensure that destructor of c has no effect
    }

    ~Field() { deallocBlock_(); }

    void *getBlockPtr() override { return static_cast<void *>(block_); }

    const void *getBlockPtr() const override
    {
        return static_cast<const void *>(block_);
    }

    size_t getBlockBytes() const override { return bytes_; }

    /// @brief Copy assignment operator
    Field &operator=(const Field &c)
    {
        if (this != &c) {
            copyBlock_(c.block_); // deep copy
        }
        return *this;
    }

    /// @brief Move assignment operator
    Field &operator=(Field &&c)
    {
        if (this != &c) {
            deallocBlock_();
            copyBlockShallow_(c.block_);
            c.setNull_(); // ensure that destructor of c has no effect
        }
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

protected:
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
    /// @param src (pointer to first block data element)
    void copyBlock_(const void *src)
    {
        if (block_ != src) {
            std::memcpy(block_, src, bytes_);
        }
    }

    /// @brief Shallow copy of an external source (required for move semantics
    ///        and field proxies)
    ///
    /// @param src (pointer to first block data element)
    void copyBlockShallow_(void *src)
    {
        block_ = static_cast<DataType *>(src);
        bytes_ = blk_alloc_.getBytes(1);
    }

    /// @brief Clear block memory bitwise
    void clearBlock_()
    {
        if (block_ != nullptr) {
            std::memset(block_, 0, bytes_);
        }
    }
};

/// @brief Field proxy type (never owns block memory).  Arithmetic operations
///        and assignment operate on the underlying data.
///
/// @tparam TField (underlying field type)
template <typename TField>
class FieldProxy : public TField
{
public:
    using FieldType = TField;

    FieldProxy() = delete;

    /// @brief Base constructor
    ///
    /// @param alloc (switch for block memory allocation)
    FieldProxy(void *block_ptr) : TField(false)
    {
        assert(block_ptr != nullptr &&
               "A field proxy can only be constructed with a valid address");
        this->block_ = static_cast<typename TField::DataType *>(block_ptr);
        this->bytes_ = this->blk_alloc_.getBytes(1);
    }

    /// @brief Copy constructor for a field proxy
    FieldProxy(const FieldProxy &c) : TField(false)
    {
        this->copyBlockShallow_(c.block_);
    }

    /// @brief Copy constructor for underlying field type
    FieldProxy(const TField &c) : TField(false)
    {
        this->copyBlockShallow_(const_cast<void *>(c.getBlockPtr()));
    }

    /// @brief Move semantics are not permitted for a field proxy
    FieldProxy(FieldProxy &&c) = delete;

    /// @brief Move semantics are not permitted for a field proxy
    FieldProxy(TField &&c) = delete;

    ~FieldProxy() { this->setNull_(); }

    /// @brief Copy assignment operator for field proxy (deep copy)
    FieldProxy &operator=(const FieldProxy &c)
    {
        if (this != &c) {
            this->copyBlock_(c.block_);
        }
        return *this;
    }

    /// @brief Copy assignment operator for underlying field type
    FieldProxy &operator=(const TField &c)
    {
        if (this != &c) {
            this->copyBlockShallow_(const_cast<void *>(c.getBlockPtr()));
        }
        return *this;
    }

    /// @brief Move semantics are not permitted for a field proxy
    FieldProxy &operator=(FieldProxy &&c) = delete;

    /// @brief Move semantics are not permitted for a field proxy
    FieldProxy &operator=(TField &&c) = delete;

    ////////////////////////////////////////////////////////////////////////////
    // TODO: [fabianw@mavt.ethz.ch; 2019-04-13] arithmetic block operators
    // (outside of class)
    ////////////////////////////////////////////////////////////////////////////
};

/// @brief Cell centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          typename BlockAlloc = AlignedBlockAllocator<DataType, BDX, BDY, BDZ>>
using FieldCell = Field<BlockAlloc, DataMapping::Cell, Cubism::Dir::Any>;

/// @brief Node centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          typename BlockAlloc =
              AlignedBlockAllocator<DataType, BDX + 1, BDY + 1, BDZ + 1>>
using FieldNode = Field<BlockAlloc, DataMapping::Node, Cubism::Dir::Any>;

/// @brief X-Face centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          typename BlockAlloc =
              AlignedBlockAllocator<DataType, BDX + 1, BDY, BDZ>>
using FieldFaceX = Field<BlockAlloc, DataMapping::Face, Cubism::Dir::X>;

/// @brief Y-Face centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          typename BlockAlloc =
              AlignedBlockAllocator<DataType, BDX, BDY + 1, BDZ>>
using FieldFaceY = Field<BlockAlloc, DataMapping::Face, Cubism::Dir::Y>;

/// @brief Z-Face centered data field type for single block
///
/// @tparam DataType
/// @tparam BDX (Block cell dimension X)
/// @tparam BDY (Block cell dimension Y)
/// @tparam BDZ (Block cell dimension Z)
/// @tparam BlockAlloc (block allocator type)
template <typename DataType,
          size_t BDX,
          size_t BDY,
          size_t BDZ,
          typename BlockAlloc =
              AlignedBlockAllocator<DataType, BDX, BDY, BDZ + 1>>
using FieldFaceZ = Field<BlockAlloc, DataMapping::Face, Cubism::Dir::Z>;

// XXX: [fabianw@mavt.ethz.ch; 2019-04-17] How to treat duplicate nodes
// and faces for DataMapping::Node or DataMapping::Face? (Treat them differently
// in Labs is one possibility, means there is duplicate data.  Problematic for
// reduction operations)

template <typename TField>
class FieldCartesian
{
public:
    FieldCartesian(size_t nbx, size_t nby = 1, size_t nbz = 1)
        : block_(nullptr), bytes_(0), nblocks_({nbx, nby, nbz})
    {
        const size_t nblocks = nblocks_[0] * nblocks_[1] * nblocks_[2];
        assert(nblocks > 0);
        allocBlocks_(nblocks);
        initProxies_(nblocks);
    }

    // TODO: [fabianw@mavt.ethz.ch; 2019-04-17] missing constructors

    virtual ~FieldCartesian() { deallocBlock_(); }

    using DataType = typename TField::DataType;
    using FieldType = TField;
    using ProxyType = FieldProxy<TField>;

    size_t getBytes() const { return bytes_; }
    size_t getNBlocks() const
    {
        return nblocks_[0] * nblocks_[1] * nblocks_[2];
    }
    size_t getNBlocksX() const { return nblocks_[0]; }
    size_t getNBlocksY() const { return nblocks_[1]; }
    size_t getNBlocksZ() const { return nblocks_[2]; }
    const std::array<size_t, 3> &getNBlocksArray() const { return nblocks_; }

    std::vector<ProxyType> &getBlocks() { return proxies_; }
    const std::vector<ProxyType> &getBlocks() const { return proxies_; }

    ProxyType &getBlock(size_t ix, size_t iy = 0, size_t iz = 0)
    {
        assert(ix < nblocks_[0]);
        assert(iy < nblocks_[1]);
        assert(iz < nblocks_[2]);
        return proxies_[ix + nblocks_[0] * (iy + nblocks_[1] * iz)];
    }
    const ProxyType &getBlock(size_t ix, size_t iy = 0, size_t iz = 0) const
    {
        assert(ix < nblocks_[0]);
        assert(iy < nblocks_[1]);
        assert(iz < nblocks_[2]);
        return proxies_[ix + nblocks_[0] * (iy + nblocks_[1] * iz)];
    }

    DataType *getData() { return block_; }
    const DataType *getData() const { return block_; }

private:
    using AllocType = typename TField::AllocType;

    DataType *block_;     // pointer to first element of first block
    size_t bytes_;        // number of bytes pointed to by block_
    AllocType blk_alloc_; // block allocator

    // Cartesian topology of blocks
    const std::array<size_t, 3> nblocks_;

    // Proxy field container
    std::vector<ProxyType> proxies_;

    /// @brief Allocate memory for multiple blocks
    ///
    /// @return True if success
    bool allocBlocks_(size_t nblocks)
    {
        block_ = blk_alloc_.allocate(nblocks);
        bytes_ = blk_alloc_.getBytes(nblocks);
        assert(block_ != nullptr && "Block allocator returned NULL address");
        if (block_ == nullptr) {
            return false;
        }
        return true;
    }

    /// @brief Deallocate all blocks
    void deallocBlock_()
    {
        if (block_ != nullptr) {
            blk_alloc_.deallocate(block_);
            bytes_ = 0;
        }
    }

    /// @brief Create a list block field proxies for interface with block data.
    ///
    /// @param nblocks
    void initProxies_(size_t nblocks)
    {
        const ptrdiff_t block_size =
            AllocType::BlockDimX * AllocType::BlockDimY * AllocType::BlockDimZ;

        proxies_.reserve(nblocks);
        for (size_t b = 0; b < nblocks; ++b) {
            void *bptr = static_cast<void *>(block_ + b * block_size);
            ProxyType fp(bptr);
            proxies_.push_back(fp);
        }
    }
};

NAMESPACE_END(BlockField)
NAMESPACE_END(Cubism)

#endif /* BLOCKFIELD_H_6LD9GMJA */
