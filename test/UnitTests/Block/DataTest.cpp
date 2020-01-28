// File       : DataTest.cpp
// Created    : Mon Dec 30 2019 06:03:31 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data test for various allocator types
// Copyright 2019 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/Data.h"
#include "Cubism/Alloc/AlignedBlockAllocator.h"
#include "Cubism/Core/Index.h"
#include "gtest/gtest.h"
#include <utility>

namespace
{
using namespace Cubism;

template <typename T, template <typename> class TAlloc, size_t DIM>
using CellData = Block::Data<T, EntityType::Cell, DIM, TAlloc<T>>;
template <typename T, template <typename> class TAlloc, size_t DIM>
using NodeData = Block::Data<T, EntityType::Node, DIM, TAlloc<T>>;

template <typename T, template <typename> class TAlloc, size_t DIM>
void runTest()
{
    using IRange = Core::IndexRange<DIM>;
    using MIndex = typename IRange::MultiIndex;

    using CData = CellData<T, TAlloc, DIM>;
    using NData = NodeData<T, TAlloc, DIM>;

    MIndex cells(16);
    MIndex nodes(16 + 1);
    IRange cell_domain(cells);
    IRange node_domain(nodes);

    CData cdata(cell_domain);
    NData ndata(node_domain);

    EXPECT_TRUE(cdata.getIndexRange() == cell_domain);
    EXPECT_EQ(cdata.getDataElementBytes(), sizeof(T));

    EXPECT_EQ(cdata.getBlockSize(), cells.prod());
    EXPECT_EQ(ndata.getBlockSize(), nodes.prod());

    // data view
    CData cdata_view(cdata, CData::MemoryOwner::No);
    EXPECT_FALSE(cdata_view.isMemoryOwner());
    EXPECT_TRUE(cdata.isMemoryOwner());
    EXPECT_EQ(cdata.getBlockPtr(), cdata_view.getBlockPtr());

    // copy construction
    NData ndata_copy(ndata);
    CData cdata_view1(cdata_view);
    EXPECT_NE(ndata.getBlockPtr(), ndata_copy.getBlockPtr());
    EXPECT_EQ(cdata_view1.getBlockPtr(), cdata_view.getBlockPtr());

    // copy assignment
    CData cdata1(cell_domain);
    cdata1 = cdata;           // deep copy
    cdata_view = cdata_view1; // shallow copy
    EXPECT_NE(cdata.getBlockPtr(), cdata1.getBlockPtr());
    EXPECT_EQ(cdata_view1.getBlockPtr(), cdata_view.getBlockPtr());

    // move construction
    CData cdata_move(std::move(cdata));
    EXPECT_EQ(cdata_view.getBlockPtr(), cdata_move.getBlockPtr());
    EXPECT_NE(cdata_move.getBlockPtr(), cdata.getBlockPtr());

    // move assignment
    NData ndata_view(ndata, NData::MemoryOwner::No);
    NData ndata_move(ndata); // deep copy
    EXPECT_NE(ndata_view.getBlockPtr(), ndata_move.getBlockPtr());
    ndata_move = std::move(ndata);
    EXPECT_EQ(ndata_view.getBlockPtr(), ndata_move.getBlockPtr());
    EXPECT_NE(ndata_view.getBlockPtr(), ndata.getBlockPtr());

    // data access
    CData ref(cdata1, CData::MemoryOwner::No);
    for (size_t i = 0; i < ref.getBlockSize(); ++i) {
        ref[i] = static_cast<T>(i);
    }

    CData test(cdata_move, CData::MemoryOwner::No);
    const auto r = test.getIndexRange();
    if (2 == DIM) { // test operator() for DIM == 2
        for (size_t j = 0; j < r.sizeDim(1); ++j) {
            for (size_t i = 0; i < r.sizeDim(0); ++i) {
                const auto val = static_cast<T>(i + r.sizeDim(0) * j);
                test(i, j) = val;
                EXPECT_EQ(test(i, j), val); // const access
            }
        }
    } else {
        for (size_t i = 0; i < test.getBlockSize(); ++i) {
            const MIndex p = r.getMultiIndex(i);
            test[p] = static_cast<T>(i);
        }
    }
    if (DIM > 3) {
        EXPECT_THROW(
            {
                try {
                    test(0);
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("Data::operator(): You can not call this "
                                 "method for DIM > 3",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    test(0) = 0;
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("Data::operator(): You can not call this "
                                 "method for DIM > 3",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
    }

    for (size_t i = 0; i < ref.getBlockSize(); ++i) {
        const MIndex p = r.getMultiIndex(i);
        EXPECT_EQ(ref[i], test[p]);
    }
}

TEST(Data, AlignedBlockAllocator)
{
    runTest<float, AlignedBlockAllocator, 1>();
    runTest<double, AlignedBlockAllocator, 2>();
    runTest<int, AlignedBlockAllocator, 4>();
}
} // namespace
