// File       : FieldLabTest.cpp
// Created    : Wed Feb 12 2020 05:13:50 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data lab test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldLab.h"
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
    using FieldLab = Block::FieldLab<Field>;
    using Stencil = typename FieldLab::StencilType;

    MIndex elements(16);
    IRange element_domain(elements);
    Field f(element_domain);
    T k = 0;
    for (auto &c : f) {
        c = k;
        k += 1;
    }
    auto fields = [&](const MIndex &) -> Field & { return f; };

    FieldLab flab;
    const Stencil s(-2, 3, TENSORIAL);
    flab.allocate(s, f.getIndexRange());
    flab.loadData(MIndex(0), fields); // periodic

    EXPECT_EQ(flab.getActiveStencil().getBegin(), s.getBegin());
    EXPECT_EQ(flab.getActiveStencil().getEnd(), s.getEnd());
    EXPECT_EQ(flab.getActiveStencil().isTensorial(), s.isTensorial());
    EXPECT_EQ(flab.getActiveStencil().isTensorial(), TENSORIAL);
    EXPECT_EQ(flab.getActiveRange().getBegin(), f.getIndexRange().getBegin());
    EXPECT_EQ(flab.getActiveRange().getEnd(), f.getIndexRange().getEnd());
    const MIndex n_max = f.getIndexRange().getExtent();
    const MIndex s_begin = s.getBegin();
    const MIndex s_end = s.getEnd();
    using Index = typename MIndex::DataType;
    for (const auto &p : flab) {
        EXPECT_EQ(flab[p], f[p]); // inner data
    }

    if (!TENSORIAL) { // ghosts
        for (const auto &p : flab) {
            if (0 == p.prod() || 0 == (p - n_max + MIndex(1)).prod()) {
                for (size_t i = 0; i < DIM; ++i) {
                    if (0 == p[i]) { // left periodic boundary
                        for (Index si = s_begin[i]; si < 0; ++si) {
                            const MIndex sm = p + si * MIndex::getUnitVector(i);
                            const MIndex sp =
                                p + (si + n_max[i]) * MIndex::getUnitVector(i);
                            EXPECT_EQ(flab[sm], flab[sp]);
                        }
                    }
                    if (0 == p[i] - (n_max[i] - 1)) { // right periodic boundary
                        for (Index si = 1; si < s_end[i]; ++si) {
                            const MIndex sm =
                                p + (si - n_max[i]) * MIndex::getUnitVector(i);
                            const MIndex sp = p + si * MIndex::getUnitVector(i);
                            EXPECT_EQ(flab[sm], flab[sp]);
                        }
                    }
                }
            }
        }
    } else {
        const IRange inner_range = flab.getActiveRange();
        const MIndex inner_extent = inner_range.getExtent();
        const MIndex ghost_begin = flab.getActiveLabRange().getBegin();
        for (const auto &p : flab.getActiveLabRange()) {
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
            EXPECT_EQ(flab[q], flab[s]);
        }
    }
}

TEST(FieldLab, Ghosts)
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

TEST(FieldLab, Reuse)
{
    using Field = Block::Field<int, Cubism::EntityType::Cell, 2>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using FieldLab = Block::FieldLab<Field>;
    using Stencil = typename FieldLab::StencilType;

    MIndex elements(16);
    IRange element_domain(elements);
    Field f(element_domain);

    FieldLab flab;
    const Stencil ssmall(-1, 2);
    const Stencil slarge(-3, 4);
    const Stencil shuge(-128, 129);
    flab.allocate(slarge, f.getIndexRange());
    const void *p0 = flab.getBlockPtr();
    const size_t b0 = flab.getBlockBytes();
    const size_t a0 = reinterpret_cast<size_t>(p0) % CUBISM_ALIGNMENT;
    EXPECT_EQ(a0, 0);
    flab.allocate(ssmall, f.getIndexRange());
    const void *p1 = flab.getBlockPtr();
    const size_t b1 = flab.getBlockBytes();
    const size_t a1 = reinterpret_cast<size_t>(p1) % CUBISM_ALIGNMENT;
    EXPECT_EQ(a1, 0);
    EXPECT_EQ(p0, p1);
    EXPECT_EQ(b0, b1);
    flab.allocate(shuge, f.getIndexRange());
    const void *p2 = flab.getBlockPtr();
    const size_t b2 = flab.getBlockBytes();
    const size_t a2 = reinterpret_cast<size_t>(p2) % CUBISM_ALIGNMENT;
    EXPECT_EQ(a2, 0);
    EXPECT_NE(p2, p1);
    EXPECT_NE(b2, b1);
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREFERENCED_FORMAL_VARIABLE

template <size_t Rank,
          Cubism::EntityType Entity,
          size_t Dim,
          size_t Comp = 0,
          size_t FDir = 0>
void testIndexer()
{
    using IRange = typename Core::IndexRange<Dim>;
    using MIndex = typename IRange::MultiIndex;

    struct MyFieldState {
        MIndex index;
    };
    using Field = typename Block::
        FieldTypeFactory<char, Rank, Entity, Dim, MyFieldState>::Type;
    using FContainer = Block::FieldContainer<Field>;

    MIndex elements(8);
    IRange element_domain(elements);
    IRange block_range(MIndex(3)); // 3 ^ Dim blocks
    FContainer fields;
    for (auto i : block_range) {
        Field *f = new Field(element_domain);
        f->getState().index = i;
        fields.pushBack(f);
    }

    using Indexer = Block::PeriodicIndexFunctor<FContainer, Field::Class, Rank>;
    Indexer i2f(fields, block_range, Comp, FDir);
    IRange block_mirror(MIndex(9));
    const MIndex block_extent = block_range.getExtent();
    const MIndex shift(-3);
    for (auto i : block_mirror) {
        const MIndex q = i + shift;
        MIndex s(q);
        for (auto &v : s) {
            v = (v + 3) % 3; // 3 as in block_range(3) above
        }
        EXPECT_EQ(i2f(q).getState().index, s);
    }
}

TEST(FieldLab, PeriodicIndexer)
{
    testIndexer<0, EntityType::Cell, 3>();
    testIndexer<0, EntityType::Node, 3>();
    testIndexer<0, EntityType::Face, 3, 0, 0>();
    testIndexer<0, EntityType::Face, 3, 0, 1>();
    testIndexer<0, EntityType::Face, 3, 0, 2>();
    testIndexer<1, EntityType::Cell, 3, 0>();
    testIndexer<1, EntityType::Cell, 3, 1>();
    testIndexer<1, EntityType::Cell, 3, 2>();
    testIndexer<1, EntityType::Node, 3, 0>();
    testIndexer<1, EntityType::Node, 3, 1>();
    testIndexer<1, EntityType::Node, 3, 2>();
    testIndexer<1, EntityType::Face, 3, 0, 0>();
    testIndexer<1, EntityType::Face, 3, 0, 1>();
    testIndexer<1, EntityType::Face, 3, 0, 2>();
    testIndexer<1, EntityType::Face, 3, 1, 0>();
    testIndexer<1, EntityType::Face, 3, 1, 1>();
    testIndexer<1, EntityType::Face, 3, 1, 2>();
    testIndexer<1, EntityType::Face, 3, 2, 0>();
    testIndexer<1, EntityType::Face, 3, 2, 1>();
    testIndexer<1, EntityType::Face, 3, 2, 2>();
}

TEST(FieldLab, Interface)
{
    using Field = Block::Field<int, Cubism::EntityType::Cell, 2>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using FieldLab = Block::FieldLab<Field>;
    using Stencil = typename FieldLab::StencilType;
    using DataType = typename Field::DataType;

    MIndex elements(16);
    IRange element_domain(elements);
    Field f(element_domain);
    DataType k = 0;
    for (auto &c : f) {
        c = k;
        k += 1;
    }

    FieldLab flab;
    const auto bytes0 = flab.getMemoryFootprint();
    EXPECT_EQ(bytes0.allocated, 0);
    EXPECT_EQ(bytes0.used, 0);

    try {
        auto &f1 = flab.getActiveField();
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ("FieldLab: no field loaded", e.what());
    }
    EXPECT_THROW(auto &f1 = flab.getActiveField(), std::runtime_error);

    try {
        const auto &f2 = flab.getActiveField();
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ("FieldLab: no field loaded", e.what());
    }
    EXPECT_THROW(const auto &f2 = flab.getActiveField(), std::runtime_error);

    const Stencil stencil(-3, 4);
    auto fields = [&](const MIndex &) -> Field & { return f; };
    flab.allocate(stencil, f.getIndexRange());
    flab.loadData(MIndex(0), fields);
    EXPECT_EQ(*flab.getInnerData(), f[0]);
    *flab.getInnerData() = 1;
    EXPECT_NE(*flab.getInnerData(), f[0]);

    EXPECT_EQ(flab.getActiveStencil().getBegin(), stencil.getBegin());
    EXPECT_EQ(flab.getActiveStencil().getEnd(), stencil.getEnd());
    EXPECT_EQ(flab.getActiveStencil().isTensorial(), stencil.isTensorial());

    const IRange arange = flab.getActiveRange();
    EXPECT_EQ(arange.getBegin(), f.getIndexRange().getBegin());
    EXPECT_EQ(arange.getEnd(), f.getIndexRange().getEnd());
    EXPECT_EQ(arange.getExtent(), f.getIndexRange().getExtent());

    auto &f3 = flab.getActiveField();
    EXPECT_EQ(&f, &f3);

    const auto &f4 = flab.getActiveField();
    EXPECT_EQ(&f, &f4);

    const IRange lrange = flab.getActiveLabRange();
    EXPECT_EQ(lrange.getBegin(), stencil.getBegin());
    EXPECT_EQ(lrange.getEnd(),
              f.getIndexRange().getExtent() + stencil.getEnd() - MIndex(1));

    const auto abytes = flab.getMemoryFootprint();
    EXPECT_EQ(abytes.allocated, flab.getBlockBytes());
    EXPECT_EQ(abytes.used, lrange.size() * sizeof(DataType));
}

DISABLE_WARNING_POP

} // namespace
