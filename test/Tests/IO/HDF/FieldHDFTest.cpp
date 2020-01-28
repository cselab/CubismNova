// File       : FieldHDFTest.cpp
// Created    : Sun Jan 26 2020 12:29:41 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Block field HDF IO
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/IO/FieldHDF.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

TEST(IO, FieldWriteHDF)
{
    // using Field = Block::CellField<float>;
    using Field = Block::NodeField<int>;
    // using Field = Block::FaceField<int>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;

    const PointType start(-1);
    const PointType end(1);
    const MIndex cells{6, 7, 8};
    Mesh m(start, end, cells, MeshIntegrity::FullMesh);
    // Field cf(m.getIndexRange(Cubism::EntityType::Cell));
    Field cf(m.getIndexRange(Cubism::EntityType::Node));
    // Field cf(m.getIndexRange(Cubism::EntityType::Face));
    int k = 0;
    for (auto &c : cf) {
        c = k++;
    }

    // all of the field
    IO::FieldWriteHDF<typename Field::DataType>("field", "field", cf, m);

    // sub-mesh selections
    const auto block = m.getSubMesh(PointType(-0.5), PointType(0.5));
    const auto slice = m.getSubMesh(PointType(-0.5), PointType{-0.5, 0.5, 0.5});
    const auto xxl = m.getSubMesh(PointType(-10.0), PointType(10.0));
    const auto lefty = m.getSubMesh(PointType{-10.0, -1.0, -10.0},
                                    PointType{10.0, -1.0, 10.0});
    const auto righty =
        m.getSubMesh(PointType{-10.0, 1.0, -10.0}, PointType{10.0, 1.0, 10.0});
    const auto linez =
        m.getSubMesh(PointType{0.0, 0.0, -10.0}, PointType{0.0, 0.0, 10.0});
    IO::FieldWriteHDF<typename Field::DataType>("block", "apples", cf, *block);
    IO::FieldWriteHDF<typename Field::DataType>("slice", "bananas", cf, *slice);
    IO::FieldWriteHDF<typename Field::DataType>("xxl", "peaches", cf, *xxl);
    IO::FieldWriteHDF<typename Field::DataType>("lefty", "oranges", cf, *lefty);
    IO::FieldWriteHDF<typename Field::DataType>("righty", "kiwis", cf, *righty);
    IO::FieldWriteHDF<typename Field::DataType>("linez", "beers", cf, *linez);
}

TEST(IO, FieldWriteReadBackHDF)
{
    using Field = Block::CellField<int>;
    // using Field = Block::NodeField<int>;
    // using Field = Block::FaceField<int>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;

    const PointType start(-1);
    const PointType end(1);
    const MIndex cells{6, 7, 8};
    Mesh m(start, end, cells, MeshIntegrity::FullMesh);
    Field src(m.getIndexRange(Cubism::EntityType::Cell));
    // Field cf(m.getIndexRange(Cubism::EntityType::Node));
    // Field cf(m.getIndexRange(Cubism::EntityType::Face));
    int k = 0;
    for (auto &c : src) {
        c = k++;
    }

    // write the field
    IO::FieldWriteHDF<typename Field::DataType>("field", "tacos", src, m);

    // read back tacos
    Field dst(m.getIndexRange(Cubism::EntityType::Cell));
    IO::FieldReadHDF<typename Field::DataType>("field", dst, m);

    k = 0;
    for (const auto &c : dst) {
        EXPECT_EQ(c, k++);
    }
}
} // namespace
