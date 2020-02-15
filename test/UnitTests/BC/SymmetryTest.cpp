// File       : SymmetryTest.cpp
// Created    : Sat Feb 15 2020 04:04:51 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Symmetry boundary test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/BC/Symmetry.h"
#include "BC/FieldAndLab.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

template <size_t dir, typename FAL, typename IRange, typename MIndex>
void check(const FAL &fal,
           const IRange &range,
           MIndex start,
           const MIndex end,
           const typename FAL::DataType s0,
           const typename FAL::DataType s1)
{
    using Index = typename MIndex::DataType;
    const auto &lab = fal.getLab();
    MIndex Nslab(range.getExtent());
    Nslab[dir] = 3;
    const IRange halo_slab(Nslab);
    Index off = -1;
    start[dir] = -3;
    for (const auto &p : halo_slab) { // side = 0
        const MIndex q = p + start;
        MIndex r(q);
        r[dir] = off - r[dir];
        EXPECT_EQ(lab[q], s0 * lab[r]);
    }
    off = 2 * end[dir] - 1;
    start[dir] = end[dir];
    for (const auto &p : halo_slab) { // side = 1
        const MIndex q = p + start;
        MIndex r(q);
        r[dir] = off - r[dir];
        EXPECT_EQ(lab[q], s1 * lab[r]);
    }
}

TEST(BC, Symmetry)
{
    using FAL = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using MIndex = typename FAL::MIndex;
    using BCVector = typename FAL::BCVector;
    using BC = BC::Symmetry<typename FAL::DataLab>;

    BCVector bcv;
    bcv.push_back(new BC(0, 0, 1));
    bcv.push_back(new BC(0, 1, -1));
    bcv.push_back(new BC(1, 0, 1));
    bcv.push_back(new BC(1, 1, -1));
    bcv.push_back(new BC(2, 0, 1));
    bcv.push_back(new BC(2, 1, -1));

    FAL fal;
    fal.loadData(&bcv);
    const auto &lab = fal.getLab();
    const auto &range = lab.getActiveRange();
    check<0>(fal, range, MIndex(0), range.getExtent(), 1, -1);
    check<1>(fal, range, MIndex(0), range.getExtent(), 1, -1);
    check<2>(fal, range, MIndex(0), range.getExtent(), 1, -1);

    for (auto bc : bcv) {
        delete bc;
    }
}

template <size_t dir, typename FAL, typename BC>
void testTensorial()
{
    using BCVector = typename FAL::BCVector;

    BCVector bcv;
    bcv.push_back(new BC(dir, 0));
    bcv.push_back(new BC(dir, 1));

    FAL fal(FAL::Tensorial::On); // tensorial stencil
    fal.loadData(&bcv);
    const auto &lab = fal.getLab();
    const auto &range = lab.getActiveLabRange();
    const auto &begin = lab.getActiveStencil().getBegin();
    check<dir>(fal, range, begin, lab.getActiveRange().getExtent(), 1, 1);

    for (auto bc : bcv) {
        delete bc;
    }
}

TEST(BC, SymmetryTensorial)
{
    using FALC = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using BCC = BC::Symmetry<typename FALC::DataLab>;
    testTensorial<0, FALC, BCC>();
    testTensorial<1, FALC, BCC>();
    testTensorial<2, FALC, BCC>();

    using FALN = FieldAndLab<int, Cubism::EntityType::Node, 3>;
    using BCN = BC::Symmetry<typename FALN::DataLab>;
    testTensorial<0, FALN, BCN>();
    testTensorial<1, FALN, BCN>();
    testTensorial<2, FALN, BCN>();
}
} // namespace
