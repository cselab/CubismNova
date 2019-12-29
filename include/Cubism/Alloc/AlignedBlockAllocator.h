// File       : AlignedBlockAllocator.h
// Created    : Sat Apr 13 2019 05:44:18 PM (+0200)
// Author     : Fabian Wermelinger
// Description: POSIX aligned block allocator
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM
#define ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM

#include "Core/Common.h"

#include <cassert>
#include <cstdlib>

NAMESPACE_BEGIN(Cubism)

/// @brief Simple aligned memory block allocator
///
/// @tparam T Data type of single block element
/// @tparam ALIGNAT Align at byte boundary
template <typename T, size_t ALIGNAT = 32>
class AlignedBlockAllocator
{
public:
    using DataType = T;

    /// @brief Allocate block memory
    ///
    /// @param nblocks Number of blocks
    ///
    /// @return Pointer to first block element
    DataType *allocate(const size_t bytes) const
    {
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
};

NAMESPACE_END(Cubism)

#endif /* ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM */
