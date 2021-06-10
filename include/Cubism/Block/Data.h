// File       : Data.h
// Created    : Sun Dec 29 2019 01:26:35 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef DATA_H_W5KVJG9U
#define DATA_H_W5KVJG9U

#include "Cubism/Alloc/AlignedBlockAllocator.h"
#include "Cubism/Common.h"
#include "Cubism/Core/Index.h"
#include <cassert>
#include <cstring>
#include <stdexcept>
#include <type_traits>

NAMESPACE_BEGIN(Cubism)
/**
 * @addtogroup Block
 * @{ */
/** @brief Namespace for block data types */
NAMESPACE_BEGIN(Block)

/**
 * @brief Type-less base class for a data block
 *
 * For byte-level operations on the data only.
 */
class DataBase
{
public:
    DataBase() = default;
    virtual ~DataBase() = default;

    DataBase(const DataBase &c) = delete;
    DataBase(DataBase &&c) = delete;

    DataBase &operator=(const DataBase &c) = delete;
    DataBase &operator=(DataBase &&c) = delete;

    /**
     * @brief Mutable pointer to address of first data element
     * @return Address of first data block byte
     */
    virtual void *getBlockPtr() { return nullptr; }

    /**
     * @brief Immutable pointer to address of first data element
     * @return Address of first data block byte
     */
    virtual const void *getBlockPtr() const { return nullptr; }

    /**
     * @brief Number of bytes occupied by block
     * @return Number of block bytes
     */
    virtual size_t getBlockBytes() const { return 0; }

    /**
     * @brief Number of bytes occupied by a single data element
     * @return Number of block bytes
     */
    virtual size_t getDataElementBytes() const { return 0; }

    /**
     * @brief Number of data elements contained in block
     * @return Number of block data elements
     */
    virtual size_t getBlockSize() const { return 0; }
};

/**
 * @brief Generic block data
 * @tparam T Data type
 * @tparam Entity Entity type
 * @tparam DIM Data dimensionality
 * @tparam BlockAlloc Allocator type
 *
 * @rst
 * Data type that manages memory allocation and data access for the specified
 * index range spanned by the data.
 *
 * .. note::
 *    The memory allocation for a block may be larger than the  minimum required
 *    data specified by the index range.
 *
 * @endrst
 */
template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename BlockAlloc = AlignedBlockAllocator<T>>
class Data : public DataBase
{
public:
    using BaseType = DataBase;
    using AllocType = BlockAlloc;
    using IndexRangeType = Core::IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using Index = typename MultiIndex::DataType;
    using DataType = typename BlockAlloc::DataType;
    static_assert(std::is_same<DataType, T>::value,
                  "Block allocator data type does not match type T");

    /** @brief Specifies memory ownership */
    enum class MemoryOwner { No = 0, Yes };

    /** @brief Byte utilization compound */
    struct BlockBytes {
        size_t allocated; // allocated bytes
        size_t used;      // used bytes
    };

    /** @brief Specifies entity type of data */
    static constexpr Cubism::EntityType EntityType = Entity;

    /**
     * @brief Base constructor for single block allocation
     * @param r Index range of data
     *
     * The index range defines the spatial dimensionality of the data.  The
     * block memory is touched by the thread that calls the constructor.
     */
    explicit Data(const IndexRangeType &r)
        : BaseType(), range_(r), owner_(MemoryOwner::Yes),
          external_memory_(false), block_(nullptr), bytes_(0)
    {
        if (static_cast<bool>(owner_)) {
            allocBlock_();
            clearBlock_(); // NUMA touch
        }
    }

    /**
     * @brief General purpose copy-constructor mainly for data views.
     * @param c Data block to copy from
     * @param owner Memory ownership
     *
     * The block memory is touched by the thread that calls the constructor.
     */
    Data(const Data &c, const MemoryOwner owner)
        : BaseType(), range_(c.range_), owner_(owner), external_memory_(false),
          block_(nullptr), bytes_(0)
    {
        if (static_cast<bool>(owner_)) {
            allocBlock_();
            copyBlockDeep_(c); // NUMA touched
        } else {
            bytes_ = c.bytes_; // explicitly set here
            copyBlockShallow_(c);
        }
    }

    /**
     * @brief Low-level constructor for externally managed memory
     * @param r Index range of data pointed to by ptr
     * @param ptr Block data pointer to first element
     * @param bytes Number of bytes of block data.
     *
     * The block owns the memory but does not deallocate it at destruction.
     */
    Data(const IndexRangeType &r, DataType *ptr, const size_t bytes)
        : BaseType(), range_(r), owner_(MemoryOwner::Yes),
          external_memory_(true), block_(ptr), bytes_(bytes)
    {
    }

    /**
     * @brief Copy constructor
     * @param c Block data to copy from
     *
     * Depending on ownership of the memory, copy is deep or shallow. Relevant
     * for views.
     */
    Data(const Data &c)
        : BaseType(), range_(c.range_), owner_(c.owner_),
          external_memory_(false), block_(nullptr), bytes_(0)
    {
        if (static_cast<bool>(owner_)) {
            allocBlock_();
            copyBlockDeep_(c); // NUMA touched
        } else {
            bytes_ = c.bytes_; // explicitly set here
            copyBlockShallow_(c);
        }
    }

    /**
     * @brief Move constructor
     * @param c Block data to copy from
     *
     * Move constructions are always shallow.  Ownership of the memory is
     * inherited from the source.
     */
    Data(Data &&c) noexcept
        : BaseType(), range_(std::move(c.range_)), owner_(c.owner_),
          external_memory_(c.external_memory_), block_(nullptr), bytes_(0)
    {
        bytes_ = c.bytes_; // explicitly set here
        copyBlockShallow_(c);
        c.setNull_(); // ensure that destructor of c has no effect
    }

    /**
     * @brief Virtual destructor
     */
    ~Data() override
    {
        if (external_memory_) {
            // if memory allocation is managed externally, nothing will be
            // deallocated in this destructor
            this->setNull_();
        }
        if (static_cast<bool>(owner_)) {
            deallocBlock_();
        }
    }

    // Virtual interface of base class
    void *getBlockPtr() override { return static_cast<void *>(block_); }
    const void *getBlockPtr() const override
    {
        return static_cast<const void *>(block_);
    }
    size_t getBlockBytes() const override { return bytes_; }
    size_t getDataElementBytes() const override { return sizeof(DataType); }
    size_t getBlockSize() const override { return range_.size(); }

    /**
     * @brief Get pointer to data
     * @return Pointer to first data element
     */
    DataType *getData() { return block_; }

    /**
     * @brief Get pointer to data
     * @return ``const`` pointer to first data element
     */
    const DataType *getData() const { return block_; }

    /**
     * @brief Copy assignment operator
     * @param c Block data to assign from
     * @return This instance with updated data
     *
     * Depending on ownership of memory, copy is deep or shallow. Relevant for
     * views.
     */
    Data &operator=(const Data &c)
    {
        if (this != &c) {
            if (static_cast<bool>(owner_)) {
                copyBlockDeep_(c);
            } else {
                copyBlockShallow_(c);
            }
        }
        return *this;
    }

    /**
     * @brief Move assignment operator
     * @param c
     * @return This instance with moved data
     *
     * Move constructions are always shallow.  Ownership of the memory is
     * inherited from the source.
     */
    Data &operator=(Data &&c)
    {
        if (this != &c) {
            if (external_memory_) {
                this->setNull_();
            }
            if (static_cast<bool>(owner_)) {
                deallocBlock_();
                bytes_ = c.bytes_; // explicitly set here
            }
            copyBlockShallow_(c);
            c.setNull_(); // ensure that destructor of c has no effect
        }
        return *this;
    }

    /**
     * @brief Linear data access
     * @param i Local flat index
     * @return Reference to data element
     */
    DataType &operator[](size_t i)
    {
        assert(i < this->getBlockSize() && "Linear index out of bounds");
        return block_[i];
    }

    /**
     * @brief Linear data access
     * @param i Local flat index
     * @return ``const`` reference to data element
     */
    const DataType &operator[](size_t i) const
    {
        assert(i < this->getBlockSize() && "Linear index out of bounds");
        return block_[i];
    }

    /**
     * @brief Linear data access
     * @param p Local multi-dimensional index
     * @return Reference to data element
     */
    DataType &operator[](const MultiIndex &p)
    {
        assert(range_.isIndex(p));
        if (1 == IndexRangeType::Dim) {
            return block_[p[0]];
        } else if (2 == IndexRangeType::Dim) {
            return block_[p[0] + range_.sizeDim(0) * p[1]];
        } else if (3 == IndexRangeType::Dim) {
            return block_[p[0] + range_.sizeDim(0) *
                                     (p[1] + range_.sizeDim(1) * p[2])];
        } else {
            return this->operator[](range_.getFlatIndex(p));
        }
    }

    /**
     * @brief Linear data access
     * @param p Local multi-dimensional index
     * @return ``const`` reference to data element
     */
    const DataType &operator[](const MultiIndex &p) const
    {
        assert(range_.isIndex(p));
        if (1 == IndexRangeType::Dim) {
            return block_[p[0]];
        } else if (2 == IndexRangeType::Dim) {
            return block_[p[0] + range_.sizeDim(0) * p[1]];
        } else if (3 == IndexRangeType::Dim) {
            return block_[p[0] + range_.sizeDim(0) *
                                     (p[1] + range_.sizeDim(1) * p[2])];
        } else {
            return this->operator[](range_.getFlatIndex(p));
        }
    }

    /**
     * @brief Classic data access
     * @param ix Index for first dimension
     * @param iy Index for second dimension
     * @param iz Index for third dimension
     * @return Reference to data element
     *
     * This operator is only supported for dimensions 1, 2 and 3.  Accessing
     * data with this operator has less latency than the multi-index sibling.
     */
    DataType &operator()(const Index ix, const Index iy = 0, const Index iz = 0)
    {
        assert((ix >= 0) && (ix < static_cast<Index>(range_.sizeDim(0))));
        if (1 == IndexRangeType::Dim) {
            return block_[ix];
        } else if (2 == IndexRangeType::Dim) {
            assert((iy >= 0) && (iy < static_cast<Index>(range_.sizeDim(1))));
            return block_[ix + range_.sizeDim(0) * iy];
        } else if (3 == IndexRangeType::Dim) {
            assert((iy >= 0) && (iy < static_cast<Index>(range_.sizeDim(1))));
            assert((iz >= 0) && (iz < static_cast<Index>(range_.sizeDim(2))));
            return block_[ix +
                          range_.sizeDim(0) * (iy + range_.sizeDim(1) * iz)];
        }
        throw std::runtime_error("FieldLab: operator() not supported");
    }

    /**
     * @brief Classic data access
     * @param ix Index for first dimension
     * @param iy Index for second dimension
     * @param iz Index for third dimension
     * @return ``const`` reference to data element
     *
     * This operator is only supported for dimensions 1, 2 and 3.  Accessing
     * data with this operator has less latency than the multi-index sibling.
     */
    const DataType &
    operator()(const Index ix, const Index iy = 0, const Index iz = 0) const
    {
        assert((ix >= 0) && (ix < static_cast<Index>(range_.sizeDim(0))));
        if (1 == IndexRangeType::Dim) {
            return block_[ix];
        } else if (2 == IndexRangeType::Dim) {
            assert((iy >= 0) && (iy < static_cast<Index>(range_.sizeDim(1))));
            return block_[ix + range_.sizeDim(0) * iy];
        } else if (3 == IndexRangeType::Dim) {
            assert((iy >= 0) && (iy < static_cast<Index>(range_.sizeDim(1))));
            assert((iz >= 0) && (iz < static_cast<Index>(range_.sizeDim(2))));
            return block_[ix +
                          range_.sizeDim(0) * (iy + range_.sizeDim(1) * iz)];
        }
        throw std::runtime_error("FieldLab: operator() not supported");
    }

    /**
     * @brief Deep copy data
     * @param c Block data to copy from
     *
     * Copies data unconditionally.
     */
    void copyData(const Data &c) { copyBlockDeep_(c); }

    /**
     * @brief Get index range
     * @return Index range spanned by the data
     */
    IndexRangeType getIndexRange(const size_t = 0) const { return range_; }

    /**
     * @brief Test for memory ownership
     * @return True if memory is owned by this instance
     */
    bool isMemoryOwner() const { return static_cast<bool>(owner_); }

    /**
     * @brief Get memory ownership
     * @return Enumeration type describing the memory ownership
     */
    MemoryOwner getMemoryOwnership() const { return owner_; }

    /**
     * @brief Get byte utilization of block
     * @return Structure of byte usage for this instance
     */
    virtual BlockBytes getMemoryFootprint() const
    {
        BlockBytes bb = {};
        bb.allocated = bytes_;
        bb.used = range_.size() * sizeof(DataType);
        return bb;
    }

protected:
    IndexRangeType range_;       // range of DIM-dimensional data
    const MemoryOwner owner_;    // owns memory
    const bool external_memory_; // true if memory allocation is external
    DataType *block_;            // pointer to first block element
    size_t bytes_;               // number of bytes pointed to by block_
    BlockAlloc blk_alloc_;       // block allocator

    /**
     * @brief Allocate single block memory
     * @return True if success
     */
    virtual bool allocBlock_()
    {
        bytes_ = range_.size() * sizeof(DataType);
        block_ = blk_alloc_.allocate(bytes_);
        assert(block_ != nullptr &&
               "Invalid memory address for block allocation");
        if (block_ == nullptr) {
            return false;
        }
        return true;
    }

    /**
     * @brief Deallocate single block memory
     */
    void deallocBlock_()
    {
        if (block_ != nullptr) {
            blk_alloc_.deallocate(block_);
            block_ = nullptr;
            bytes_ = 0;
        }
    }

    /**
     * @brief Set all addresses to NULL
     */
    void setNull_()
    {
        block_ = nullptr;
        bytes_ = 0;
    }

    /**
     * @brief Deep copy of an external source
     * @param rhs Data block to be copied
     */
    void copyBlockDeep_(const Data &rhs)
    {
        static_assert(std::is_pod<DataType>::value, "DataType is not POD");
        assert(bytes_ == rhs.bytes_ &&
               "copyBlockDeep: Number of block bytes are not identical");
        if (block_ != rhs.block_) {
            std::memcpy(block_, rhs.block_, bytes_);
        }
    }

    /**
     * @brief Shallow copy of an external source
     * @param rhs Data block to be copied
     */
    void copyBlockShallow_(const Data &rhs)
    {
        assert(bytes_ == rhs.bytes_ &&
               "copyBlockShallow: Number of block bytes are not identical");
        block_ = rhs.block_;
        bytes_ = rhs.bytes_;
    }

    /**
     * @brief Clear block memory bitwise
     */
    void clearBlock_()
    {
        static_assert(std::is_pod<DataType>::value, "DataType is not POD");
        if (block_ != nullptr) {
            std::memset(block_, 0, bytes_);
        }
    }
};

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename BlockAlloc>
constexpr Cubism::EntityType Data<T, Entity, DIM, BlockAlloc>::EntityType;

NAMESPACE_END(Block)
/**  @} */
NAMESPACE_END(Cubism)

#endif /* DATA_H_W5KVJG9U */
