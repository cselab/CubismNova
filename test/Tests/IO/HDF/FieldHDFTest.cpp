// File       : FieldHDFTest.cpp
// Created    : Sun Jan 26 2020 12:29:41 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Block field HDF IO
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "IO/FieldHDF.h"
#include "Block/Field.h"
#include "Core/Index.h"
#include "Mesh/StructuredUniform.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

TEST(IO, FieldWriteHDF)
{
    using CellField = Block::CellField<float>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;

    const PointType start(-1);
    const PointType end(1);
    const MIndex cells{6, 7, 8};
    Mesh m(start, end, cells, MeshIntegrity::FullMesh);
    CellField cf(m.getIndexRange(Cubism::EntityType::Cell));
    int k = 0;
    for (auto &c : cf) {
        c = k++;
    }
    IO::FieldWriteHDF<typename CellField::DataType>("test", "test", cf, m);
}
} // namespace
