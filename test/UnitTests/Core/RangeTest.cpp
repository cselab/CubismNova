// File       : RangeTest.cpp
// Created    : Sun Dec 29 2019 09:31:34 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Generic ranges test
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Core/Range.h"
#include "gtest/gtest.h"

namespace
{
// some test types
using Range = Cubism::Core::Range<float, 3>;
using Point = typename Range::PointType;

// Construction of Ranges
TEST(Range, Construction)
{
    EXPECT_EQ(Range::Dim, 3);

    Range r0(2);
    Range r1(0, 2);

    Point p0(0);
    Point p1(2);
    Range r2(p1);
    Range r3(p0, p1);

    EXPECT_EQ(r0, r1);
    EXPECT_EQ(r1, r2);
    EXPECT_EQ(r2, r3);

    Range r4(r3);
    Range r5(3);
    r5 = r4;
    EXPECT_EQ(r3, r4);
    EXPECT_EQ(r4, r5);

    Range r6(Range(2));
    r5 = Range(2);
    EXPECT_EQ(r4, r5);
    EXPECT_EQ(r5, r6);

    EXPECT_THROW(
        {
            try {
                Range r7(p1, p0);
            } catch (const std::runtime_error &e) {
                EXPECT_STREQ(
                    "RangeConstruction: begin_ must be smaller than end_",
                    e.what());
                throw;
            }
        },
        std::runtime_error);
}

// Set and get of Ranges
TEST(Range, SetGet)
{
    Point p0(1);
    Point p1(2);
    Point p2(3);
    EXPECT_THROW(
        {
            try {
                Range r0(p0);
                r0.setBegin(p1);
            } catch (const std::runtime_error &e) {
                EXPECT_STREQ("RangeSetBegin: begin_ must be smaller than end_",
                             e.what());
                throw;
            }
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            try {
                Range r0(p1, p2);
                r0.setEnd(p0);
            } catch (const std::runtime_error &e) {
                EXPECT_STREQ("RangeSetEnd: begin_ must be smaller than end_",
                             e.what());
                throw;
            }
        },
        std::runtime_error);

    Range r0(p0);
    Range r1(p1, p2);
    r0.setEnd(p2);
    r0.setBegin(p1);
    EXPECT_EQ(r0, r1);
    EXPECT_EQ(r0.getBegin(), p1);
    EXPECT_EQ(r0.getEnd(), p2);

    EXPECT_EQ(r0.getExtent(), p2 - p1);
    EXPECT_EQ(r0.getVolume(), (p2 - p1).prod());
}

// Utilities
TEST(Range, Utils)
{
    Point p0(1);
    Point p1(2);
    Point p2(3);
    Point p3(4);
    Range r0(p3);
    Range r1(p1, p2);
    Range r2(p2);
    Range r3(p1, p3);
    Range r4(p0);

    EXPECT_TRUE(r0.isContained(r1));
    EXPECT_TRUE(r0.isContained(p0));
    EXPECT_TRUE(r2.isIntersecting(r1));
    EXPECT_TRUE(r3.isIntersecting(r2));
    EXPECT_FALSE(r4.isIntersecting(r1));
    EXPECT_FALSE(r1.isIntersecting(r4));
}

} // namespace
