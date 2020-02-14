// File       : BaseTest.cpp
// Created    : Fri Feb 14 2020 07:05:48 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Boundary conditions base class test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/BC/Base.h"
#include "BC/FieldAndLab.h"
#include "gtest/gtest.h"

namespace
{
using namespace Cubism;

TEST(BC, Base)
{
    using FAL = FieldAndLab<int, Cubism::EntityType::Cell, 3>;
    using BCVector = typename FAL::BCVector;
    using BC = BC::Base<typename FAL::DataLab>;

    BCVector bcv;
    bcv.push_back(new BC(0, 0));
    bcv.push_back(new BC(0, 1));
    bcv.push_back(new BC(1, 0));
    bcv.push_back(new BC(1, 1));
    bcv.push_back(new BC(2, 0));
    bcv.push_back(new BC(2, 1));

    FAL fal;
    fal.loadData(&bcv);

    using MIndex = typename FAL::MIndex;
    using Index = typename MIndex::DataType;
    const auto &lab = fal.getLab();
    const MIndex n_max = lab.getActiveRange().getExtent();
    const MIndex s_begin = lab.getActiveStencil().getBegin();
    const MIndex s_end = lab.getActiveStencil().getEnd();
    for (const auto &p : lab) {
        if (0 == p.prod() || 0 == (p - n_max + MIndex(1)).prod()) {
            for (size_t i = 0; i < 3; ++i) {
                if (0 == p[i]) { // left periodic boundary
                    for (Index si = s_begin[i]; si < 0; ++si) {
                        const MIndex sm = p + si * MIndex::getUnitVector(i);
                        const MIndex sp =
                            p + (si + n_max[i]) * MIndex::getUnitVector(i);
                        EXPECT_EQ(lab[sm], lab[sp]);
                    }
                }
                if (0 == p[i] - (n_max[i] - 1)) { // right periodic boundary
                    for (Index si = 1; si < s_end[i]; ++si) {
                        const MIndex sm =
                            p + (si - n_max[i]) * MIndex::getUnitVector(i);
                        const MIndex sp = p + si * MIndex::getUnitVector(i);
                        EXPECT_EQ(lab[sm], lab[sp]);
                    }
                }
            }
        }
    }

    for (auto bc : bcv) {
        delete bc;
    }
}
} // namespace
