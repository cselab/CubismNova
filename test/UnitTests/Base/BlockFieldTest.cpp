// File       : BlockFieldTest.cpp
// Created    : Fri May 03 2019 04:44:17 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Unit tests for Cubism/Base/BlockField.h
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Base/BlockField.h"
#include "Core/Common.h"
#include "gtest/gtest.h"

namespace
{
// block size
#define BSX 32
#define BSY 64
#define BSZ 128

using DataType = double;
using FieldCell = Cubism::BlockField::FieldCell<DataType, BSX, BSY, BSZ>;
using FieldNode = Cubism::BlockField::FieldNode<DataType, BSX, BSY, BSZ>;
using FieldFaceX = Cubism::BlockField::FieldFaceX<DataType, BSX, BSY, BSZ>;
using FieldFaceY = Cubism::BlockField::FieldFaceY<DataType, BSX, BSY, BSZ>;
using FieldFaceZ = Cubism::BlockField::FieldFaceZ<DataType, BSX, BSY, BSZ>;

template <typename TField>
void setFieldValue(TField &f, const typename TField::DataType v)
{
    for (size_t i = 0; i < f.getBlockSize(); ++i) {
        f[i] = v;
    }
}

template <typename TField>
typename TField::DataType sumField(const TField &f)
{
    typename TField::DataType sum = 0;
    for (size_t i = 0; i < f.getBlockSize(); ++i) {
        sum += f[i];
    }
    return sum;
}

// Static members
TEST(FieldTest, Static)
{
    ASSERT_EQ(FieldCell::BlockDimX, BSX);
    ASSERT_EQ(FieldCell::BlockDimY, BSY);
    ASSERT_EQ(FieldCell::BlockDimZ, BSZ);
    ASSERT_EQ(FieldCell::MapClass, Cubism::DataMapping::Cell);
    ASSERT_EQ(FieldCell::Dir, Cubism::Dir::Any);

    ASSERT_EQ(FieldNode::BlockDimX, BSX + 1);
    ASSERT_EQ(FieldNode::BlockDimY, BSY + 1);
    ASSERT_EQ(FieldNode::BlockDimZ, BSZ + 1);
    ASSERT_EQ(FieldNode::MapClass, Cubism::DataMapping::Node);
    ASSERT_EQ(FieldNode::Dir, Cubism::Dir::Any);

    ASSERT_EQ(FieldFaceX::BlockDimX, BSX + 1);
    ASSERT_EQ(FieldFaceX::BlockDimY, BSY);
    ASSERT_EQ(FieldFaceX::BlockDimZ, BSZ);
    ASSERT_EQ(FieldFaceX::MapClass, Cubism::DataMapping::Face);
    ASSERT_EQ(FieldFaceX::Dir, Cubism::Dir::X);

    ASSERT_EQ(FieldFaceY::BlockDimX, BSX);
    ASSERT_EQ(FieldFaceY::BlockDimY, BSY + 1);
    ASSERT_EQ(FieldFaceY::BlockDimZ, BSZ);
    ASSERT_EQ(FieldFaceY::MapClass, Cubism::DataMapping::Face);
    ASSERT_EQ(FieldFaceY::Dir, Cubism::Dir::Y);

    ASSERT_EQ(FieldFaceZ::BlockDimX, BSX);
    ASSERT_EQ(FieldFaceZ::BlockDimY, BSY);
    ASSERT_EQ(FieldFaceZ::BlockDimZ, BSZ + 1);
    ASSERT_EQ(FieldFaceZ::MapClass, Cubism::DataMapping::Face);
    ASSERT_EQ(FieldFaceZ::Dir, Cubism::Dir::Z);
}

TEST(FieldTest, ConstructorAndAssignment)
{
    FieldCell fc0; // default
    EXPECT_NE(fc0.getBlockPtr(), nullptr);
    EXPECT_EQ(fc0.getBlockBytes(), BSX * BSY * BSZ * sizeof(DataType));

    FieldCell fc1(false); // no memory allocated
    EXPECT_EQ(fc1.getBlockPtr(), nullptr);
    EXPECT_EQ(fc1.getBlockBytes(), 0);

    FieldCell fc2(fc0); // copy
    EXPECT_NE(fc2.getBlockPtr(), fc0.getBlockPtr());
    EXPECT_EQ(fc2.getBlockBytes(), fc0.getBlockBytes());

    const void *const pfc0 = fc0.getBlockPtr();
    FieldCell fc3(std::move(fc0)); // move construct
    EXPECT_EQ(fc3.getBlockPtr(), pfc0);
    EXPECT_EQ(fc0.getBlockPtr(), nullptr);
    EXPECT_EQ(fc0.getBlockBytes(), 0);

    setFieldValue(fc3, 1);
    fc2 = fc3; // assign
    EXPECT_EQ(sumField(fc2), fc2.getBlockSize());

    const void *const pfc2 = fc2.getBlockPtr();
    fc3 = std::move(fc2); // move assign
    EXPECT_EQ(fc3.getBlockPtr(), pfc2);
    EXPECT_EQ(fc2.getBlockPtr(), nullptr);
    EXPECT_EQ(fc2.getBlockBytes(), 0);
}

TEST(FieldTest, BasePointer)
{
    FieldCell fc;
    const typename FieldCell::BaseType *pc = &fc;
    EXPECT_EQ(pc->getBlockSize(), fc.getBlockSize());
}

TEST(ProxyFieldTest, ConstructorAndAssignment)
{
    // TODO: [fabianw@mavt.ethz.ch; 2019-09-20] revise
    using FieldProxy = Cubism::BlockField::FieldProxy<FieldCell>;

    FieldCell fc0; // memory carrier

    FieldProxy fp0(fc0); // construct from field
    EXPECT_EQ(fp0.getBlockPtr(), fc0.getBlockPtr());
    EXPECT_EQ(fp0.getBlockBytes(), fc0.getBlockBytes());

    // FieldProxy fp1(fc0.getBlockPtr()); // construct with block pointer
    // EXPECT_EQ(fp1.getBlockPtr(), fc0.getBlockPtr());
    // EXPECT_EQ(fp1.getBlockBytes(), fc0.getBlockBytes());

    FieldProxy fp2(fp0); // copy construct from other proxy
    EXPECT_EQ(fp2.getBlockPtr(), fp0.getBlockPtr());
    EXPECT_EQ(fp2.getBlockBytes(), fp0.getBlockBytes());

    FieldCell fc1; // another carrier
    setFieldValue(fc1, 1);
    FieldProxy fp3(fc1);
    fp2 = fp3; // deep copy field pointed to by fp3 into field pointed to by fp2
    EXPECT_EQ(sumField(fp2), fp2.getBlockSize());

    fp2 = fc1; // shallow copy field fc1 into proxy fp1
    EXPECT_EQ(sumField(fp2), fp2.getBlockSize());
}

} // namespace
