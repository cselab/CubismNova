// File       : VersionTest.cpp
// Created    : Thu Jan 16 2020 01:42:50 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cubism version string tests
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Util/Version.h"
#include "gtest/gtest.h"
#include <iostream>

namespace
{
using namespace Cubism;

TEST(Version, Stdout)
{
    std::cout << "CubismVersion     = " << Util::CubismVersion << std::endl;
}
} // namespace
