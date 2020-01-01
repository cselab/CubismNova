// File       : FieldTest.cpp
// Created    : Mon Dec 30 2019 11:35:34 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block field test
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Block/Field.h"
#include "Core/Index.h"
#include "gtest/gtest.h"

#include "Alloc/AlignedBlockAllocator.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace
{
using namespace Cubism;

template <typename T, template <typename> class TAlloc, size_t DIM>
using CellData = Block::Data<T, DataMapping::Cell, DIM, TAlloc<T>>;
template <typename T, template <typename> class TAlloc, size_t DIM>
using NodeData = Block::Data<T, DataMapping::Node, DIM, TAlloc<T>>;
template <typename T, template <typename> class TAlloc, size_t DIM>
using FaceData = Block::Data<T, DataMapping::Face, DIM, TAlloc<T>>;

TEST(Field, Construction)
{
    using CellField = Block::Field<CellData<int, AlignedBlockAllocator, 3>>;
    using FieldState = typename CellField::FieldStateType;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using FieldView = Block::FieldView<CellField>;

    MIndex cells(16);
    IRange cell_domain(cells);

    { // default
        CellField cf(cell_domain);
        EXPECT_TRUE(cf.isScalar());
        EXPECT_EQ(cf.getRank(), 0);
        EXPECT_EQ(cf.getComp(), 0);

        FieldState fs;
        fs.rank = 1;
        fs.comp = 1;
        CellField cf1(cell_domain, fs);
        EXPECT_FALSE(cf1.isScalar());
        EXPECT_EQ(cf1.getRank(), 1);
        EXPECT_EQ(cf1.getComp(), 1);
    }

    { // low-level
        CellField cf(cell_domain);
        CellField cf1(cf, CellField::BaseType::MemoryOwner::Yes);
        FieldState fs;
        fs.rank = cf.getRank();
        fs.comp = cf.getComp();
        CellField cf2(cell_domain, cf.getData(), cf.getBlockBytes(), &fs);

        EXPECT_EQ(cf.isMemoryOwner(), cf1.isMemoryOwner());
        EXPECT_NE(cf.getBlockPtr(), cf1.getBlockPtr());
        EXPECT_EQ(cf.isScalar(), cf1.isScalar());
        EXPECT_EQ(cf.getRank(), cf1.getRank());
        EXPECT_EQ(cf.getComp(), cf1.getComp());

        EXPECT_NE(cf.isMemoryOwner(), cf2.isMemoryOwner());
        EXPECT_EQ(cf.getBlockPtr(), cf2.getBlockPtr());
        EXPECT_EQ(cf.isScalar(), cf2.isScalar());
        EXPECT_EQ(cf.getRank(), cf2.getRank());
        EXPECT_EQ(cf.getComp(), cf2.getComp());
    }

    { // copy construction
        CellField cf(cell_domain);
        CellField cf_copy(cf);
        FieldView cf_view(cf);

        EXPECT_EQ(cf.isMemoryOwner(), cf_copy.isMemoryOwner());
        EXPECT_NE(cf.getBlockPtr(), cf_copy.getBlockPtr());
        EXPECT_EQ(cf.isScalar(), cf_copy.isScalar());
        EXPECT_EQ(cf.getRank(), cf_copy.getRank());
        EXPECT_EQ(cf.getComp(), cf_copy.getComp());

        EXPECT_NE(cf.isMemoryOwner(), cf_view.isMemoryOwner());
        EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
        EXPECT_EQ(cf.isScalar(), cf_view.isScalar());
        EXPECT_EQ(cf.getRank(), cf_view.getRank());
        EXPECT_EQ(cf.getComp(), cf_view.getComp());
    }

    { // copy assignment
        CellField cf(cell_domain);

        FieldState fs;
        fs.rank = 1;
        fs.comp = 1;
        CellField cf1(cell_domain, fs);
        CellField cf_copy(cf);
        FieldView cf_view(cf);

        EXPECT_EQ(cf.getRank(), 0);
        EXPECT_EQ(cf.getComp(), 0);
        EXPECT_EQ(cf_copy.getRank(), 0);
        EXPECT_EQ(cf_copy.getComp(), 0);
        EXPECT_EQ(cf_view.getRank(), 0);
        EXPECT_EQ(cf_view.getComp(), 0);
        EXPECT_EQ(cf1.getRank(), 1);
        EXPECT_EQ(cf1.getComp(), 1);
        EXPECT_NE(cf.getBlockPtr(), cf1.getBlockPtr());
        EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
        EXPECT_EQ(&cf.getState(), &cf_view.getState());

        cf = cf1; // deep copy
        EXPECT_EQ(cf.getRank(), 1);
        EXPECT_EQ(cf.getComp(), 1);
        EXPECT_EQ(cf_copy.getRank(), 0);
        EXPECT_EQ(cf_copy.getComp(), 0);
        EXPECT_EQ(cf_view.getRank(), 1);
        EXPECT_EQ(cf_view.getComp(), 1);
        EXPECT_NE(cf.getBlockPtr(), cf1.getBlockPtr());
        EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
        EXPECT_EQ(&cf.getState(), &cf_view.getState());
        EXPECT_NE(&cf_copy.getState(), &cf_view.getState());
    }

    { // move construction
        CellField cf(cell_domain);
        FieldView cfv(cf);

        EXPECT_TRUE(cf.isScalar());
        EXPECT_EQ(cf.getRank(), 0);
        EXPECT_EQ(cf.getComp(), 0);

        CellField cf_move(std::move(cf));
        EXPECT_EQ(cf.getBlockPtr(), nullptr);
        EXPECT_EQ(&cf.getState(), nullptr);
        EXPECT_TRUE(cf_move.isScalar());
        EXPECT_EQ(cf_move.getRank(), 0);
        EXPECT_EQ(cf_move.getComp(), 0);
        EXPECT_EQ(cfv.getRank(), 0);
        EXPECT_EQ(cfv.getComp(), 0);
        EXPECT_EQ(cfv.getBlockPtr(), cf_move.getBlockPtr());
        EXPECT_EQ(&cfv.getState(), &cf_move.getState());
    }

    { // move assignment
        CellField cf(cell_domain);
        FieldView cfv(cf);

        EXPECT_TRUE(cf.isScalar());
        EXPECT_EQ(cf.getRank(), 0);
        EXPECT_EQ(cf.getComp(), 0);

        CellField cf_move(cell_domain);
        cf_move = std::move(cf);
        EXPECT_EQ(cf.getBlockPtr(), nullptr);
        EXPECT_EQ(&cf.getState(), nullptr);
        EXPECT_TRUE(cf_move.isScalar());
        EXPECT_EQ(cf_move.getRank(), 0);
        EXPECT_EQ(cf_move.getComp(), 0);
        EXPECT_EQ(cfv.getRank(), 0);
        EXPECT_EQ(cfv.getComp(), 0);
        EXPECT_EQ(cfv.getBlockPtr(), cf_move.getBlockPtr());
        EXPECT_EQ(&cfv.getState(), &cf_move.getState());
    }
}

TEST(Field, Interface)
{
    using FaceField = Block::Field<FaceData<double, AlignedBlockAllocator, 3>>;
    using FieldState = typename FaceField::FieldStateType;
    using IRange = typename FaceField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex faces(16);
    IRange face_domain(faces);
    FaceField ff(face_domain);
    const FieldState fs = ff.getState();

    EXPECT_TRUE(ff.isScalar());
    EXPECT_EQ(ff.getRank(), 0);
    EXPECT_EQ(ff.getComp(), 0);
    EXPECT_EQ(ff.getRank(), fs.rank);
    EXPECT_EQ(ff.getComp(), fs.comp);
}

TEST(Field, Iterator)
{
    using CellField = Block::Field<CellData<float, AlignedBlockAllocator, 4>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    CellField cf(cell_domain);

    std::fill(cf.begin(), cf.end(), 1);

    typename CellField::DataType sum = 0;
    for (const auto c : cf) {
        sum += c;
    }
    EXPECT_EQ(sum, cell_domain.size());
}

TEST(Field, View)
{
    using NodeField = Block::Field<NodeData<int, AlignedBlockAllocator, 2>>;
    using IRange = typename NodeField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    using FieldView = Block::FieldView<NodeField>;
    MIndex nodes(16);
    IRange node_domain(nodes);

    NodeField nf(node_domain);
    FieldView nfv(nf);

    EXPECT_NE(&nf, &nfv);
    EXPECT_EQ(nf.getBlockPtr(), nfv.getBlockPtr());
    EXPECT_TRUE(nf.isMemoryOwner());
    EXPECT_FALSE(nfv.isMemoryOwner());
}

// TODO: [fabianw@mavt.ethz.ch; 2020-01-01]
TEST(Field, Arithmetic) {}

} // namespace
