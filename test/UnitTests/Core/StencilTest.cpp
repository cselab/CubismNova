// File       : StencilTest.cpp
// Created    : Wed Feb 12 2020 04:54:46 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Test stencil type
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Core/Stencil.h"
#include "gtest/gtest.h"

namespace
{
template <size_t DIM>
void runTest()
{
    using Stencil = Cubism::Core::Stencil<DIM>;
    using PointType = typename Stencil::PointType;

    Stencil s0(-2, 3);
    EXPECT_EQ(s0.getBegin(), PointType(-2));
    EXPECT_EQ(s0.getEnd(), PointType(3));
    EXPECT_FALSE(s0.isTensorial());

    PointType b1(-6);
    PointType e1(3);
    Stencil s1(b1, e1, true);
    EXPECT_EQ(s1.getBegin(), b1);
    EXPECT_EQ(s1.getEnd(), e1);
    EXPECT_TRUE(s1.isTensorial());

    try {
        Stencil s2(0, 0);
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ("Stencil: end_ must be > 0", e.what());
    }
    EXPECT_THROW(Stencil s2(0, 0), std::runtime_error);

    try {
        Stencil s3(1, 4);
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ("Stencil: begin_ must be <= 0", e.what());
    }
    EXPECT_THROW(Stencil s3(1, 4), std::runtime_error);
}

TEST(Stencil, Interface)
{
    runTest<1>();
    runTest<2>();
    runTest<3>();
}
} // namespace
