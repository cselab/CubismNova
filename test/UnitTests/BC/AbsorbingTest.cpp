// File       : AbsorbingTest.cpp
// Created    : Fri Feb 14 2020 07:58:04 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Zeroth-Order absorbing boundary test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/BC/Absorbing.h"
#include "BC/FieldAndLab.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

template <size_t dir, typename FAL, typename IRange, typename MIndex>
void check(const FAL &fal, const IRange &range, MIndex start, const MIndex end)
{
    using Index = typename MIndex::DataType;
    const auto &lab = fal.getLab();
    MIndex Nslab(range.getExtent());
    Nslab[dir] = 3;
    const IRange halo_slab(Nslab);
    Index src = 0;
    start[dir] = -3;
    for (const auto &p : halo_slab) { // side = 0
        const MIndex q = p + start;
        MIndex r(q);
        r[dir] = src;
        EXPECT_EQ(lab[q], lab[r]);
    }
    src = end[dir] - 1;
    start[dir] = end[dir];
    for (const auto &p : halo_slab) { // side = 1
        const MIndex q = p + start;
        MIndex r(q);
        r[dir] = src;
        EXPECT_EQ(lab[q], lab[r]);
    }
}

TEST(BC, Absorbing)
{
    using FAL = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using MIndex = typename FAL::MIndex;
    using BCVector = typename FAL::BCVector;
    using BC = BC::Absorbing<typename FAL::DataLab>;

    BCVector bcv;
    bcv.push_back(new BC(0, 0));
    bcv.push_back(new BC(0, 1));
    bcv.push_back(new BC(1, 0));
    bcv.push_back(new BC(1, 1));
    bcv.push_back(new BC(2, 0));
    bcv.push_back(new BC(2, 1));

    FAL fal;
    fal.loadData(&bcv);
    const auto &lab = fal.getLab();
    const auto &range = lab.getActiveRange();
    check<0>(fal, range, MIndex(0), range.getExtent());
    check<1>(fal, range, MIndex(0), range.getExtent());
    check<2>(fal, range, MIndex(0), range.getExtent());

    for (auto bc : bcv) {
        delete bc;
    }
}

template <size_t dir, typename FAL, typename BC>
void testDirection()
{
    using BCVector = typename FAL::BCVector;

    BCVector bcv;
    bcv.push_back(new BC(dir, 0, true)); // tensorial boundary
    bcv.push_back(new BC(dir, 1, true)); // tensorial boundary

    FAL fal(FAL::Tensorial::On);
    fal.loadData(&bcv);
    const auto &lab = fal.getLab();
    const auto &range = lab.getActiveLabRange();
    const auto &begin = lab.getActiveStencil().getBegin();
    check<dir>(fal, range, begin, lab.getActiveRange().getExtent());

    for (auto bc : bcv) {
        delete bc;
    }
}

TEST(BC, AbsorbingTensorial)
{
    using FAL = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using BC = BC::Absorbing<typename FAL::DataLab>;
    testDirection<0, FAL, BC>();
    testDirection<1, FAL, BC>();
    testDirection<2, FAL, BC>();
}
} // namespace
