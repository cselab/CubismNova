// File       : DataHDFTest.cpp
// Created    : Wed Jun 09 2021 05:41:14 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Block data HDF IO
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/Field.h"
#include "Cubism/IO/Data.h"
#include "gtest/gtest.h"
#include <numeric>

namespace
{
using namespace Cubism;

TEST(IO, DataWriteUniformHDF)
{
    using CellField = Block::CellField<int>;
    using NodeField = Block::NodeField<int>;
    using IRange = typename CellField::IndexRangeType;

    IRange elements(8);

    { // cell field
        CellField f(elements);
        std::iota(f.begin(), f.end(), 0);
        IO::DataWriteUniformHDF<int>(
            "ucdata",
            "data",
            static_cast<typename CellField::BlockDataType>(f));
    }
    { // node field
        NodeField f(elements);
        std::iota(f.begin(), f.end(), 0);
        IO::DataWriteUniformHDF<float>(
            "undata",
            "data",
            static_cast<typename NodeField::BlockDataType>(f));
    }
}
} // namespace
