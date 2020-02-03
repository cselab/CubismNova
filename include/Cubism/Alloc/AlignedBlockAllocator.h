// File       : AlignedBlockAllocator.h
// Created    : Sat Apr 13 2019 05:44:18 PM (+0200)
// Author     : Fabian Wermelinger
// Description: POSIX aligned block allocator
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM
#define ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM

#include "Cubism/Common.h"
#include <cassert>
#include <cstdlib>

NAMESPACE_BEGIN(Cubism)

/**
 * @brief Simple aligned memory block allocator
 * @tparam T Data type of single block element
 */
template <typename T>
class AlignedBlockAllocator
{
public:
    using DataType = T;
    static constexpr size_t Alignment = CUBISM_ALIGNMENT;

    /**
     * @brief Allocate block memory
     * @param bytes Minimum number of bytes
     * @return Pointer to first block element
     *
     * The actual allocated memory may be larger than the requested number of
     * bytes.
     */
    DataType *allocate(size_t &bytes) const
    {
        void *block = nullptr;
        // ensure byte block is an integer multiple of CUBISM_ALIGNMENT
        bytes = ((bytes + CUBISM_ALIGNMENT - 1) / CUBISM_ALIGNMENT) *
                CUBISM_ALIGNMENT;
        int ret = posix_memalign(&block, CUBISM_ALIGNMENT, bytes);
        assert(ret == 0 && block != nullptr &&
               "posix_memalign returned NULL address");
        (void)ret;
        return static_cast<DataType *>(block);
    }

    /**
     * @brief Deallocate block memory
     * @param block Pointer to first block element
     */
    void deallocate(DataType *block) const
    {
        if (block != nullptr) {
            free(block);
        }
    }
};

template <typename T>
constexpr size_t AlignedBlockAllocator<T>::Alignment;

NAMESPACE_END(Cubism)

#endif /* ALIGNEDBLOCKALLOCATOR_H_ARTMHJKM */
