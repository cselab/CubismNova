// File       : INIParserTest.cpp
// Created    : Tue Dec 24 2019 12:21:53 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Parse test for .ini parser
// Copyright 2019 ETH Zurich. All Rights Reserved.

#include "Cubism/Util/INIParser.h"
#include "gtest/gtest.h"
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
TEST(INIParser, FileRead)
{
    Cubism::Util::INIParser p("INIFiles/main.ini");
    ASSERT_FALSE(p.parseError());

    for (const auto e : p.fileErrors()) {
        ASSERT_EQ(e.second, 0);
    }
}

TEST(INIParser, Interface)
{
    Cubism::Util::INIParser p("INIFiles/main.ini");
    std::cout << p;

    // defined sections
    EXPECT_TRUE(p.hasSection("inc_a"));
    EXPECT_TRUE(p.hasSection("inc_aaa"));
    EXPECT_TRUE(p.hasSection("include"));
    EXPECT_TRUE(p.hasSection("main"));
    EXPECT_TRUE(p.hasSection("test"));

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
    EXPECT_FALSE(p.hasSection("sparta"));      // section does not exist
    EXPECT_FALSE(p.hasValue("sparta", "300")); // section and key do not exist

    // runtime throws
    {
        EXPECT_THROW(
            {
                try {
                    p.get("sparta", "300");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ(
                        "get: key=300 in section=sparta does not exist",
                        e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    p.get("test", "noval");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("get: key=noval in section=test has no value",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    p.get("test", "empty");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("get: key=empty in section=test has no value",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);

        // scalar values
        EXPECT_NO_THROW(p.getString("test", "good"));
        EXPECT_NO_THROW(p.getString("test", "bad"));
        EXPECT_NO_THROW(p.getInteger("test", "good"));
        EXPECT_NO_THROW(p.getReal("test", "good"));
        EXPECT_NO_THROW(p.getBoolean("test", "good"));
        EXPECT_THROW(
            {
                try {
                    p.getInteger("test", "bad");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("getInteger: can not convert 'ouch' to "
                                 "integer for key=bad in section=test",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    p.getReal("test", "bad");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("getReal: can not convert 'ouch' to "
                                 "floating point for key=bad in section=test",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    p.getBoolean("test", "bad");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("getBoolean: can not convert 'ouch' to "
                                 "boolean for key=bad in section=test",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);

        // arrays
        EXPECT_NO_THROW(p.getStringArray("test", "array"));
        EXPECT_THROW(
            {
                try {
                    p.getIntegerArray("test", "array");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("getIntegerArray: can not convert 'ouch' to "
                                 "integer for key=array in section=test",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    p.getRealArray("test", "array");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("getRealArray: can not convert 'ouch' to "
                                 "floating point for key=array in section=test",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    p.getBooleanArray("test", "array");
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("getBooleanArray: can not convert '2' to "
                                 "boolean for key=array in section=test",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
    }

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
