// File       : CartesianTest.cpp
// Created    : Sun Jan 05 2020 09:27:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian Grid test
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Grid/Cartesian.h"
#include "Core/Index.h"
#include "Mesh/StructuredUniform.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

TEST(Cartesian, Construction)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;

    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;

    // using MeshHull = typename Mesh::MeshHull;
    // using PointType = typename Mesh::PointType;
    // using Entity = typename Mesh::EntityType;
    // using Range = typename Mesh::RangeType;
    using Grid = Grid::Cartesian<float, Mesh>;

    const MIndex nblocks(2);
    const MIndex block_cells(8);
    Grid g(nblocks, block_cells);
}

} // namespace
