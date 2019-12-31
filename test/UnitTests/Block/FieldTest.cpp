// File       : FieldTest.cpp
// Created    : Mon Dec 30 2019 11:35:34 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block field test
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Block/Field.h"
#include "Core/Index.h"
#include "gtest/gtest.h"

#include "Alloc/AlignedBlockAllocator.h"

#include <utility>

#include <algorithm>
#include <iostream>

namespace
{
using namespace Cubism;

template <typename T, template <typename> class TAlloc, size_t DIM>
using CellData = Block::Data<T, DataMapping::Cell, DIM, TAlloc<T>>;
template <typename T, template <typename> class TAlloc, size_t DIM>
using NodeData = Block::Data<T, DataMapping::Node, DIM, TAlloc<T>>;

TEST(Field, Iterator)
{
    using CellField3D = Block::Field<CellData<int, AlignedBlockAllocator, 3>>;
    using IRange = typename CellField3D::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    CellField3D cf(cell_domain);

    Block::FieldContainer<CellField3D> fc(3, cell_domain);

    auto &cfr = fc[Cubism::Dir::X];
    std::fill(cfr.begin(), cfr.end(), 1);

    // Block::FieldView<CellField3D> cfv(cf);
    // auto cfv(cf);

    // EXPECT_FALSE(cfv.isMemoryOwner());
    // EXPECT_EQ(cf.getBlockPtr(), cfv.getBlockPtr());

    std::fill(cf.begin(), cf.end(), 1);

    for (auto c : cf) {
        std::cout << c << std::endl;
    }
}
} // namespace
