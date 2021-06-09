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

template <Cubism::EntityType Entity>
struct Initializer {
    template <typename Field>
    void init(Field &f)
    {
        typename Field::DataType k = 0;
        for (auto &c : f) {
            c = k++;
        }
    }
};

template <>
struct Initializer<Cubism::EntityType::Face> {
    template <typename Field>
    void init(Field &f)
    {
        typename Field::DataType k = 0;
        for (auto d : f) {
            for (auto &c : *d) {
                c = k++;
            }
        }
    }
};

TEST(IO, FieldWriteHDF)
{
    using Field =
        Block::FieldTypeFactory<int, 0, Cubism::EntityType::Cell>::Type;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;

    const PointType start(-1);
    const PointType end(1);
    const MIndex cells{6, 7, 8};
    Mesh m(start, end, cells, MeshIntegrity::FullMesh);

    // sub-mesh selections
    //
    // sub-block (from -0.5 to 0.5)
    const auto blk = m.getSubMesh(PointType(-0.5), PointType(0.5));
    // interior slice at x=-0.5 (from -0.5 to 0.5)
    const auto ice = m.getSubMesh(PointType(-0.5), PointType{-0.5, 0.5, 0.5});
    // larger sub-mesh request will default to same mesh as m
    const auto xxl = m.getSubMesh(PointType(-10.0), PointType(10.0));
    // slice at y=-1 (left boundary)
    const auto lfy = m.getSubMesh(PointType{-10.0, -1.0, -10.0},
                                  PointType{10.0, -1.0, 10.0});
    // slice at y=1 (right boundary)
    const auto rty =
        m.getSubMesh(PointType{-10.0, 1.0, -10.0}, PointType{10.0, 1.0, 10.0});
    // line at (0,0) along z-direction
    const auto liz =
        m.getSubMesh(PointType{0.0, 0.0, -10.0}, PointType{0.0, 0.0, 10.0});

    { // cell field
        using DataType = typename Field::DataType;
        Field f(m.getIndexRange(Field::EntityType));
        Initializer<Field::EntityType> finit;
        finit.init(f);
        const double time = 0;
        IO::FieldWriteHDF<DataType>("call", "basket", f, m, time);
        IO::FieldWriteHDF<DataType>("cblk", "apples", f, *blk, time);
        IO::FieldWriteHDF<DataType>("cice", "bananas", f, *ice, time);
        IO::FieldWriteHDF<DataType>("cxxl", "peaches", f, *xxl, time);
        IO::FieldWriteHDF<DataType>("clfy", "oranges", f, *lfy, time);
        IO::FieldWriteHDF<DataType>("crty", "kiwis", f, *rty, time);
        IO::FieldWriteHDF<DataType>("cliz", "beers", f, *liz, time);
    }
    { // node field
        using Field =
            Block::FieldTypeFactory<int, 0, Cubism::EntityType::Node>::Type;
        using DataType = typename Field::DataType;
        Field f(m.getIndexRange(Field::EntityType));
        Initializer<Field::EntityType> finit;
        finit.init(f);
        const double time = 0;
        IO::FieldWriteHDF<DataType>("nall", "basket", f, m, time);
        IO::FieldWriteHDF<DataType>("nblk", "apples", f, *blk, time);
        IO::FieldWriteHDF<DataType>("nice", "bananas", f, *ice, time);
        IO::FieldWriteHDF<DataType>("nxxl", "peaches", f, *xxl, time);
        IO::FieldWriteHDF<DataType>("nlfy", "oranges", f, *lfy, time);
        IO::FieldWriteHDF<DataType>("nrty", "kiwis", f, *rty, time);
        IO::FieldWriteHDF<DataType>("nliz", "beers", f, *liz, time);
    }
    { // face field's
        using Field =
            Block::FieldTypeFactory<int, 0, Cubism::EntityType::Face>::Type;
        using DataType = typename Field::DataType;
        Field f(m.getIndexRange(Field::EntityType));
        Initializer<Field::EntityType> finit;
        finit.init(f);
        const double time = 0;
        for (size_t d = 0; d < Field::IndexRangeType::Dim; ++d) {
            IO::FieldWriteHDF<DataType>(
                "fall" + std::to_string(d), "basket", f, m, time, d);
            IO::FieldWriteHDF<DataType>(
                "fblk" + std::to_string(d), "apples", f, *blk, time, d);
            IO::FieldWriteHDF<DataType>(
                "fice" + std::to_string(d), "bananas", f, *ice, time, d);
            IO::FieldWriteHDF<DataType>(
                "fxxl" + std::to_string(d), "peaches", f, *xxl, time, d);
            IO::FieldWriteHDF<DataType>(
                "flfy" + std::to_string(d), "oranges", f, *lfy, time, d);
            IO::FieldWriteHDF<DataType>(
                "frty" + std::to_string(d), "kiwis", f, *rty, time, d);
            IO::FieldWriteHDF<DataType>(
                "fliz" + std::to_string(d), "beers", f, *liz, time, d);
        }
    }
}

TEST(IO, FieldWriteReadBackHDF)
{
    using Field =
        Block::FieldTypeFactory<int, 0, Cubism::EntityType::Cell>::Type;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;

    const PointType start(-1);
    const PointType end(1);
    const MIndex cells{6, 7, 8};
    Mesh m(start, end, cells, MeshIntegrity::FullMesh);
    Field src(m.getIndexRange(Field::EntityType));
    Initializer<Field::EntityType> finit;
    finit.init(src);

    // write the field
    IO::FieldWriteHDF<typename Field::DataType>("ftacos", "tacos", src, m, 0);

    // read back tacos
    Field dst(m.getIndexRange(Cubism::EntityType::Cell));
    IO::FieldReadHDF<typename Field::DataType>("ftacos", dst, m);

    int k = 0;
    for (const auto &c : dst) {
        EXPECT_EQ(c, k++);
    }
}
} // namespace
