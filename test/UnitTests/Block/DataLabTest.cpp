// File       : DataLabTest.cpp
// Created    : Wed Feb 12 2020 05:13:50 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data lab test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/DataLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

template <typename T,
          size_t DIM,
          Cubism::EntityType Entity,
          bool TENSORIAL = false>
void runTest()
{
    using Field = Block::Field<T, Entity, DIM>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using DataLab = Block::DataLab<Field>;
    using Stencil = typename DataLab::StencilType;

    MIndex elements(16);
    IRange element_domain(elements);
    Field f(element_domain);
    T k = 0;
    for (auto &c : f) {
        c = k;
        k += 1;
    }
    auto fields = [&](const MIndex &) -> const Field & { return f; };

    DataLab dlab;
    const Stencil s(-2, 3, TENSORIAL);
    dlab.allocate(s, f.getIndexRange());
    dlab.loadData(MIndex(0), fields); // periodic

    EXPECT_EQ(dlab.getActiveStencil().getBegin(), s.getBegin());
    EXPECT_EQ(dlab.getActiveStencil().getEnd(), s.getEnd());
    EXPECT_EQ(dlab.getActiveStencil().isTensorial(), s.isTensorial());
    EXPECT_EQ(dlab.getActiveStencil().isTensorial(), TENSORIAL);
    EXPECT_EQ(dlab.getActiveRange().getBegin(), f.getIndexRange().getBegin());
    EXPECT_EQ(dlab.getActiveRange().getEnd(), f.getIndexRange().getEnd());
    const MIndex n_max = f.getIndexRange().getExtent();
    const MIndex s_begin = s.getBegin();
    const MIndex s_end = s.getEnd();
    using Index = typename MIndex::DataType;
    for (const auto &p : dlab) {
        EXPECT_EQ(dlab[p], f[p]); // inner data
    }

    if (!TENSORIAL) { // ghosts
        for (const auto &p : dlab) {
            if (0 == p.prod() || 0 == (p - n_max + MIndex(1)).prod()) {
                for (size_t i = 0; i < DIM; ++i) {
                    if (0 == p[i]) { // left periodic boundary
                        for (Index si = s_begin[i]; si < 0; ++si) {
                            const MIndex sm = p + si * MIndex::getUnitVector(i);
                            const MIndex sp =
                                p + (si + n_max[i]) * MIndex::getUnitVector(i);
                            EXPECT_EQ(dlab[sm], dlab[sp]);
                        }
                    }
                    if (0 == p[i] - (n_max[i] - 1)) { // right periodic boundary
                        for (Index si = 1; si < s_end[i]; ++si) {
                            const MIndex sm =
                                p + (si - n_max[i]) * MIndex::getUnitVector(i);
                            const MIndex sp = p + si * MIndex::getUnitVector(i);
                            EXPECT_EQ(dlab[sm], dlab[sp]);
                        }
                    }
                }
            }
        }
    } else {
        const MIndex one(1);
        const IRange inner_range = dlab.getActiveRange();
        const MIndex inner_extent = inner_range.getExtent();
        const MIndex ghost_begin = dlab.getActiveLabRange().getBegin();
        for (const auto &p : dlab.getActiveLabRange()) {
            const MIndex q = p + ghost_begin;
            if (inner_range.isIndex(q)) {
                continue;
            }
            MIndex s(q); // periodic projection
            for (size_t i = 0; i < DIM; ++i) {
                if (s[i] < 0) {
                    s[i] += inner_extent[i];
                } else if (s[i] >= inner_extent[i]) {
                    s[i] = s[i] % inner_extent[i];
                }
            }
            EXPECT_EQ(dlab[q], dlab[s]);
        }
    }
}

TEST(DataLab, Ghosts)
{
    runTest<int, 1, Cubism::EntityType::Cell, false>();
    runTest<int, 1, Cubism::EntityType::Cell, true>();
    runTest<int, 2, Cubism::EntityType::Cell, false>();
    runTest<int, 2, Cubism::EntityType::Cell, true>();
    runTest<int, 3, Cubism::EntityType::Cell, false>();
    runTest<int, 3, Cubism::EntityType::Cell, true>();

    runTest<int, 2, Cubism::EntityType::Node, false>();
    runTest<int, 2, Cubism::EntityType::Node, true>();
    runTest<int, 2, Cubism::EntityType::Face, false>();
    runTest<int, 2, Cubism::EntityType::Face, true>();
}
} // namespace
