// File       : DataLabTest.cpp
// Created    : Wed Feb 12 2020 05:13:50 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data lab test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/DataLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"
#include "gtest/gtest.h"

// XXX: [fabianw@mavt.ethz.ch; 2020-02-12] debug
#include "Cubism/IO/FieldHDF.h"
#include "Cubism/Mesh/StructuredUniform.h"

namespace
{
using namespace Cubism;

template <typename T,
          size_t DIM,
          Cubism::EntityType Entity>
void runTest()
{
    using Field = Block::Field<T, Entity, DIM>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using DataLab = Block::DataLab<Field>;
    using Stencil = typename DataLab::StencilType;

    MIndex cells(16);
    IRange cell_domain(cells);
    Field f(cell_domain);
    int k = 0;
    for (auto &c : f) {
        c = k++;
    }
    auto fields = [&](const MIndex &) -> const Field & { return f; };

    DataLab dlab;
    const Stencil s(-2, 3);
    dlab.allocate(s, f.getIndexRange());
    dlab.loadData(MIndex(0), fields);

    // XXX: [fabianw@mavt.ethz.ch; 2020-02-12] debug
    using Mesh = Mesh::StructuredUniform<T, DIM>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;
    const Mesh m(PointType(1),
                 dlab.getIndexRange().getExtent(),
                 MeshIntegrity::FullMesh);
    Field flab(dlab);
    IO::FieldWriteHDF<T>("dlab", "data", flab, m, 0);
}

TEST(DataLab, Basic) { runTest<float, 2, Cubism::EntityType::Cell>(); }
} // namespace
