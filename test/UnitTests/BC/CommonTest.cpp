// File       : CommonTest.cpp
// Created    : Fri Jun 04 2021 09:31:01 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Common tests for boundary conditions
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/BC/Absorbing.h"
#include "Cubism/BC/Base.h"
#include "Cubism/BC/Dirichlet.h"
#include "Cubism/BC/Symmetry.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldLab.h"
#include "gtest/gtest.h"
#include <numeric>

namespace
{
using namespace Cubism;

template <typename Field, typename FieldLab>
void testBC(BC::Base<FieldLab> *boundary)
{
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using Stencil = typename FieldLab::StencilType;

    const Stencil s(0, 1); // empty stencil

    MIndex elements(4);
    IRange element_domain(elements);
    Field f(element_domain);
    std::iota(f.begin(), f.end(), 0);

    FieldLab flab;
    auto fields = [&](const MIndex &) -> Field & { return f; }; // periodic
    flab.allocate(s, f.getIndexRange());
    flab.loadData(MIndex(0), fields);

    boundary->operator()(flab); // apply boundary to field lab
    for (const auto &p : f.getIndexRange()) {
        EXPECT_EQ(flab[p], f[p]);
    }
}

TEST(BC, NoStencilWidth)
{
    using Field = Block::CellField<double>;
    using FieldLab = Block::FieldLab<Field>;
    {
        BC::Absorbing<FieldLab> bc(0, 0);
        testBC<Field, FieldLab>(&bc);
    }
    {
        BC::Dirichlet<FieldLab> bc(0, 0, 1.0);
        testBC<Field, FieldLab>(&bc);
    }
    {
        BC::Symmetry<FieldLab> bc(0, 0);
        testBC<Field, FieldLab>(&bc);
    }
}
} // namespace
