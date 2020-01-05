// File       : CartesianTest.cpp
// Created    : Sun Jan 05 2020 09:27:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian Grid test
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Grid/Cartesian.h"

namespace
{
using namespace Cubism;

template <typename T, template <typename> class TAlloc, size_t DIM>
using CellData = Block::Data<T, EntityType::Cell, DIM, TAlloc<T>>;

TEST(Cartesian, Construction)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;

    using CellField =
        Block::Field<CellData<double, AlignedBlockAllocator, Mesh::Dim>,
                     MyFieldState>;

    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using PointType = typename Mesh::PointType;
    using Range = typename Mesh::RangeType;

    const MIndex block_cells(8);
    const MIndex nblocks(2);
}
} // namespace
