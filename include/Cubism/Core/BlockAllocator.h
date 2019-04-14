// File       : BlockAllocator.h
// Created    : Sat Apr 13 2019 05:44:18 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Various block allocator types for use with Base/BlockField.h
//              types
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef BLOCKALLOCATOR_H_W0YCMRLB
#define BLOCKALLOCATOR_H_W0YCMRLB

#include "Core/Common.h"

#include <cassert>
#include <cstdlib>

NAMESPACE_BEGIN(Cubism)

/// @brief Simple aligned memory block allocator
///
/// @tparam T Data type of single block element
/// @tparam DIMX Number of block data elements in X dimension
/// @tparam DIMY Number of block data elements in Y dimension
/// @tparam DIMZ Number of block data elements in Z dimension
/// @tparam ALIGNAT Align at byte boundary
template <typename T,
          size_t DIMX,
          size_t DIMY,
          size_t DIMZ,
          size_t ALIGNAT = 32>
class AlignedBlockAllocator
{
public:
    using DataType = T;
    static constexpr size_t BlockDimX = DIMX;
    static constexpr size_t BlockDimY = DIMY;
    static constexpr size_t BlockDimZ = DIMZ;

    /// @brief Allocate block memory
    ///
    /// @param nblocks Number of blocks
    ///
    /// @return Pointer to first block element
    DataType *allocate(size_t nblocks) const
    {
        const size_t bytes = getBytes(nblocks);
        void *block = nullptr;
        posix_memalign(&block, ALIGNAT, bytes);
        assert(block != nullptr && "posix_memalign returned NULL address");
        return static_cast<DataType *>(block);
    }

    /// @brief Deallocate block memory
    ///
    /// @param block Reference of pointer to first block element
    void deallocate(DataType *&block) const
    {
        if (block != nullptr) {
            free(block);
            block = nullptr;
        }
    }

    /// @brief Memory footprint occupied by block(s)
    ///
    /// @param nblocks Number of blocks
    ///
    /// @return Total block memory footprint in bytes
    size_t getBytes(size_t nblocks) const
    {
        return nblocks * BlockDimX * BlockDimY * BlockDimZ * sizeof(DataType);
    }
};

NAMESPACE_END(Cubism)

#endif /* BLOCKALLOCATOR_H_W0YCMRLB */
