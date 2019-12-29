// File       : IndexTest.cpp
// Created    : Sun Dec 29 2019 09:31:02 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Test indexing space
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Core/Index.h"
#include "gtest/gtest.h"

namespace
{
using IRange = Cubism::IndexRange<3>;
using MIndex = IRange::PointType;

TEST(Index, ExtendedInterface)
{
    MIndex p0(1);
    MIndex p1(2);

    IRange r0(1);
    IRange r1(p0);
    IRange r2(2);
    IRange r3(p1);
    IRange r4(1, 2);
    IRange r5(p0, p1);

    EXPECT_EQ(r0, r1);
    EXPECT_EQ(r2, r3);
    EXPECT_EQ(r4, r5);

    EXPECT_EQ(r0.size(), r1.size());
    EXPECT_EQ(r2.size(), r3.size());
    EXPECT_EQ(r4.size(), r5.size());

    EXPECT_EQ(r4.sizeDim(0), r5.sizeDim(0));
    EXPECT_EQ(r4.sizeDim(1), r5.sizeDim(1));
    EXPECT_EQ(r4.sizeDim(2), r5.sizeDim(2));
}

} // namespace
