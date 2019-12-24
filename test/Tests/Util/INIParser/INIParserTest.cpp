// File       : INIParserTest.cpp
// Created    : Tue Dec 24 2019 12:21:53 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Parse test for .ini parser
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Util/INIParser.h"
#include "gtest/gtest.h"

#include <iostream>
#include <vector>

namespace
{
TEST(INIParser, FileRead)
{
    Cubism::Util::INIParser p("INIFiles/main.ini");
    ASSERT_FALSE(p.parseError());
    std::cout << p;

    for (const auto e : p.fileErrors()) {
        ASSERT_EQ(e.second, 0);
    }

    // defined sections
    EXPECT_TRUE(p.hasSection("include"));
    EXPECT_TRUE(p.hasSection("main"));
    EXPECT_TRUE(p.hasSection("no_value"));
    EXPECT_TRUE(p.hasSection("inc_a"));
    EXPECT_TRUE(p.hasSection("inc_aaa"));

    // check if keys have values
    EXPECT_TRUE(p.hasValue("main", "A"));
    EXPECT_TRUE(p.hasValue("main", "ivector"));
    EXPECT_TRUE(p.hasValue("main", "scalar"));
    EXPECT_TRUE(p.hasValue("main", "vector"));
    EXPECT_TRUE(p.hasValue("main", "bool"));
    EXPECT_TRUE(p.hasValue("main", "boolArray"));
    EXPECT_TRUE(p.hasValue("inc_a", "special"));
    EXPECT_TRUE(p.hasValue("inc_aaa", "special"));

    // check if has no value
    EXPECT_FALSE(p.hasValue("no_value", "dummy"));

    { // strings
        const auto str = p.getString("inc_aaa", "special");
        EXPECT_EQ(str, std::string("four more special things"));

        const auto vec = p.getStringArray("inc_aaa", "special");
        EXPECT_EQ(vec.size(), 4); // from incAaa.ini
        EXPECT_EQ(vec[0], std::string("four"));
        EXPECT_EQ(vec[1], std::string("more"));
        EXPECT_EQ(vec[2], std::string("special"));
        EXPECT_EQ(vec[3], std::string("things"));
    }

    { // integer
        const int i = p.getInteger("main", "A");
        EXPECT_EQ(i, 2); // included from 'incAaa.ini'

        const auto dvec = p.getIntegerArray("main", "vector");
        EXPECT_EQ(dvec.size(), 3); // from main.ini
        EXPECT_EQ(dvec[0], 0);
        EXPECT_EQ(dvec[1], 1);
        EXPECT_EQ(dvec[2], 2);

        const auto ivec = p.getIntegerArray("main", "ivector");
        EXPECT_EQ(ivec.size(), 4); // from main.ini
        EXPECT_EQ(ivec[0], 0);
        EXPECT_EQ(ivec[1], 1);
        EXPECT_EQ(ivec[2], 1);
        EXPECT_EQ(ivec[3], 1);
    }

    { // real
        const double d = p.getReal("main", "scalar");
        EXPECT_EQ(d, 2.0); // included from 'incB.ini'

        const float f = p.getReal("main", "scalar");
        EXPECT_EQ(f, 2.0); // included from 'incB.ini'

        const auto dvec = p.getRealArray("main", "vector");
        EXPECT_EQ(dvec.size(), 3); // from main.ini
        EXPECT_EQ(dvec[0], 0.0);
        EXPECT_EQ(dvec[1], 1.0);
        EXPECT_EQ(dvec[2], 2.012);
    }

    { // boolean
        const bool b = p.getBoolean("main", "bool");
        EXPECT_TRUE(b); // included from 'main.ini'

        const auto bvec = p.getBooleanArray("main", "boolArray");
        EXPECT_EQ(bvec.size(), 8); // from main.ini
        EXPECT_TRUE(bvec[0]);
        EXPECT_TRUE(bvec[1]);
        EXPECT_TRUE(bvec[2]);
        EXPECT_TRUE(bvec[3]);
        EXPECT_FALSE(bvec[4]);
        EXPECT_FALSE(bvec[5]);
        EXPECT_FALSE(bvec[6]);
        EXPECT_FALSE(bvec[7]);
    }
}
} // namespace
