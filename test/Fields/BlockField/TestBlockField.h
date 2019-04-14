// File       : TestBlockField.h
// Created    : Sun Apr 14 2019 03:32:49 PM (+0200)
// Author     : Fabian Wermelinger
// Description: BlockField test
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef TESTBLOCKFIELD_H_CMAHYZVB
#define TESTBLOCKFIELD_H_CMAHYZVB

#include "Base/BlockField.h"
#include "Core/Timer.h"

#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

using DataType = double;

static void tell(const std::string &msg) { std::cout << msg << '\n'; }

template <typename T>
static void setFieldValue(T &f, const typename T::DataType v)
{
    for (size_t i = 0; i < f.getBlockSize(); ++i) {
        f[i] = v;
    }
}

template <typename T>
static typename T::DataType sumField(const T &f)
{
    typename T::DataType sum = 0;
    for (size_t i = 0; i < f.getBlockSize(); ++i) {
        sum += f[i];
    }
    return sum;
}

template <typename T>
static void printStaticMember(std::ostream &s)
{
    const std::string name(typeid(T).name());
    s << name << "::BlockDimX = " << T::BlockDimX << '\n';
    s << name << "::BlockDimY = " << T::BlockDimY << '\n';
    s << name << "::BlockDimZ = " << T::BlockDimZ << '\n';
    s << name << "::MapName   = " << T::MapName << '\n';
    s << name << "::MapClass  = " << static_cast<size_t>(T::MapClass) << '\n';
    s << name << "::Dim       = " << static_cast<size_t>(T::Dim) << '\n';
}

#endif /* TESTBLOCKFIELD_H_CMAHYZVB */
