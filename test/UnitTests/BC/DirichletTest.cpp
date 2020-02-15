// File       : DirichletTest.cpp
// Created    : Sat Feb 15 2020 02:53:57 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Dirichlet boundary test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/BC/Dirichlet.h"
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
           const typename FAL::DataType v0,
           const typename FAL::DataType v1)
{
    const auto &lab = fal.getLab();
    MIndex Nslab(range.getExtent());
    Nslab[dir] = 3;
    const IRange halo_slab(Nslab);
    start[dir] = -3;
    for (const auto &p : halo_slab) { // side = 0
        EXPECT_EQ(lab[p + start], v0);
    }
    start[dir] = end[dir];
    for (const auto &p : halo_slab) { // side = 1
        EXPECT_EQ(lab[p + start], v1);
    }
}

TEST(BC, Dirichlet)
{
    using FAL = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using MIndex = typename FAL::MIndex;
    using BCVector = typename FAL::BCVector;
    using BC = BC::Dirichlet<typename FAL::DataLab>;

    BCVector bcv;
    bcv.push_back(new BC(0, 0, 10));
    bcv.push_back(new BC(0, 1, 11));
    bcv.push_back(new BC(1, 0, 12));
    bcv.push_back(new BC(1, 1, 13));
    bcv.push_back(new BC(2, 0, 14));
    bcv.push_back(new BC(2, 1, 15));

    FAL fal;
    fal.loadData(&bcv);
    const auto &lab = fal.getLab();
    const auto &range = lab.getActiveRange();
    check<0>(fal, range, MIndex(0), range.getExtent(), 10, 11);
    check<1>(fal, range, MIndex(0), range.getExtent(), 12, 13);
    check<2>(fal, range, MIndex(0), range.getExtent(), 14, 15);

    for (auto bc : bcv) {
        delete bc;
    }
}

template <size_t dir, typename FAL, typename BC>
void testTensorial(const typename FAL::DataType v0,
                   const typename FAL::DataType v1)
{
    using BCVector = typename FAL::BCVector;

    BCVector bcv;
    bcv.push_back(new BC(dir, 0, v0));
    bcv.push_back(new BC(dir, 1, v1));

    FAL fal(FAL::Tensorial::On); // tensorial stencil
    fal.loadData(&bcv);
    const auto &lab = fal.getLab();
    const auto &range = lab.getActiveLabRange();
    const auto &begin = lab.getActiveStencil().getBegin();
    check<dir>(fal, range, begin, lab.getActiveRange().getExtent(), v0, v1);

    for (auto bc : bcv) {
        delete bc;
    }
}

TEST(BC, DirichletTensorial)
{
    using FALC = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using BCC = BC::Dirichlet<typename FALC::DataLab>;
    testTensorial<0, FALC, BCC>(101, 102);
    testTensorial<1, FALC, BCC>(103, 104);
    testTensorial<2, FALC, BCC>(105, 106);

    using FALN = FieldAndLab<int, Cubism::EntityType::Node, 3>;
    using BCN = BC::Dirichlet<typename FALN::DataLab>;
    testTensorial<0, FALN, BCN>(201, 202);
    testTensorial<1, FALN, BCN>(203, 204);
    testTensorial<2, FALN, BCN>(205, 206);
}

TEST(BC, DirichletInterface)
{
    using FAL = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using BC = BC::Dirichlet<typename FAL::DataLab>;

    BC *bc = new BC(0, 0, 10);
    EXPECT_EQ(bc->getValue(), 10);
    bc->getValue() = 101;
    EXPECT_EQ(bc->getValue(), 101);

    delete bc;
}
} // namespace
