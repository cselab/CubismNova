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
template <typename T, size_t ALIGNAT = CUBISM_ALIGNMENT>
class AlignedBlockAllocator
{
public:
    using DataType = T;
    static constexpr size_t Alignment = ALIGNAT;

    /// @brief Allocate block memory
    ///
    /// @param bytes Minimum number of bytes
    ///
    /// @return Pointer to first block element
    DataType *allocate(size_t &bytes) const
    {
        void *block = nullptr;
        // ensure byte block is an integer multiple of ALIGNAT
        bytes = ((bytes + ALIGNAT - 1) / ALIGNAT) * ALIGNAT;
        posix_memalign(&block, ALIGNAT, bytes);
        assert(block != nullptr && "posix_memalign returned NULL address");
        return static_cast<DataType *>(block);
    }

    /// @brief Deallocate block memory
    ///
    /// @param block Pointer to first block element
    void deallocate(DataType *block) const
    {
        if (block != nullptr) {
            free(block);
        }
    }
};

template <typename T, size_t ALIGNAT>
constexpr size_t AlignedBlockAllocator<T, ALIGNAT>::Alignment;

NAMESPACE_END(Cubism)

#endif /* ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM */
