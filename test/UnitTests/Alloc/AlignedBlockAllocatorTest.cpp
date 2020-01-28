// File       : AlignedBlockAllocatorTest.cpp
// Created    : Sun Dec 29 2019 09:34:16 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Test aligned POSIX block allocator
// Copyright 2019 ETH Zurich. All Rights Reserved.

#include "Cubism/Alloc/AlignedBlockAllocator.h"
#include "gtest/gtest.h"

namespace
{
TEST(Alloc, AlignedBlockAllocator)
{
    using Alloc = Cubism::AlignedBlockAllocator<int>;
    using T = Alloc::DataType;

    Alloc a;
    constexpr size_t N = 10;
    size_t bytes = N * sizeof(T);
    EXPECT_FALSE(bytes % Alloc::Alignment == 0);

    T *aptr = a.allocate(bytes);

    // test alignment
    EXPECT_TRUE(bytes % Alloc::Alignment == 0);
    EXPECT_TRUE(reinterpret_cast<size_t>(aptr) % Alloc::Alignment == 0);

    // assign some values
    for (size_t i = 0; i < N; ++i) {
        aptr[i] = static_cast<int>(i);
        EXPECT_EQ(aptr[i], static_cast<int>(i));
    }

    a.deallocate(aptr);
}
} // namespace
