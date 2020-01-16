// File       : VersionTest.cpp
// Created    : Thu Jan 16 2020 01:42:50 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cubism version string tests
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Util/Version.h"
#include "gtest/gtest.h"

#include <iostream>

namespace
{
using namespace Cubism;

TEST(Version, Stdout)
{
    std::cout << "CubismVersion     = " << Util::CubismVersion << std::endl;
    std::cout << "CubismVersionHEAD = " << Util::CubismVersionHEAD << std::endl;
    std::cout << "CubismBranch      = " << Util::CubismBranch << std::endl;
}
} // namespace
