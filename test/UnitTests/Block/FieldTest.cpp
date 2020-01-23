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
#include <type_traits>
#include <utility>

namespace
{
using namespace Cubism;

template <typename T, template <typename> class TAlloc, size_t DIM>
using CellData = Block::Data<T, EntityType::Cell, DIM, TAlloc<T>>;
template <typename T, template <typename> class TAlloc, size_t DIM>
using NodeData = Block::Data<T, EntityType::Node, DIM, TAlloc<T>>;
template <typename T, template <typename> class TAlloc, size_t DIM>
using FaceData = Block::Data<T, EntityType::Face, DIM, TAlloc<T>>;

TEST(Field, Construction)
{
    using CellField = Block::Field<CellData<int, AlignedBlockAllocator, 3>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using FieldView = Block::FieldView<CellField>;

    MIndex cells(16);
    IRange cell_domain(cells);

    { // default
        CellField cf(cell_domain);
        EXPECT_TRUE(cf.isScalar());
    }

    { // copy construction
        CellField cf(cell_domain);
        CellField cf_copy(cf);
        FieldView cf_view(cf);

        EXPECT_EQ(cf.isMemoryOwner(), cf_copy.isMemoryOwner());
        EXPECT_NE(cf.getBlockPtr(), cf_copy.getBlockPtr());
        EXPECT_NE(&cf.getState(), &cf_copy.getState());
        EXPECT_EQ(cf.isScalar(), cf_copy.isScalar());

        EXPECT_NE(cf.isMemoryOwner(), cf_view.isMemoryOwner());
        EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
        EXPECT_EQ(&cf.getState(), &cf_view.getState());
        EXPECT_EQ(cf.isScalar(), cf_view.isScalar());
    }

    { // copy assignment
        CellField cf(cell_domain);
        CellField cf1(cell_domain);
        FieldView cf_view(cf);

        EXPECT_NE(cf.getBlockPtr(), cf1.getBlockPtr());
        EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
        EXPECT_EQ(&cf.getState(), &cf_view.getState());

        cf = cf1; // deep copy
        EXPECT_NE(cf.getBlockPtr(), cf1.getBlockPtr());
        EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
        EXPECT_EQ(&cf.getState(), &cf_view.getState());
    }

    { // move construction
        CellField cf(cell_domain);
        FieldView cfv(cf);

        EXPECT_TRUE(cf.isScalar());

        CellField cf_move(std::move(cf));
        EXPECT_EQ(cf.getBlockPtr(), nullptr);
        EXPECT_EQ(&cf.getState(), nullptr);
        EXPECT_TRUE(cf_move.isScalar());
        EXPECT_EQ(cfv.getBlockPtr(), cf_move.getBlockPtr());
        EXPECT_EQ(&cfv.getState(), &cf_move.getState());
    }

    { // move assignment
        CellField cf(cell_domain);
        FieldView cfv(cf);

        EXPECT_TRUE(cf.isScalar());

        CellField cf_move(cell_domain);
        cf_move = std::move(cf);
        EXPECT_EQ(cf.getBlockPtr(), nullptr);
        EXPECT_EQ(&cf.getState(), nullptr);
        EXPECT_TRUE(cf_move.isScalar());
        EXPECT_EQ(cfv.getBlockPtr(), cf_move.getBlockPtr());
        EXPECT_EQ(&cfv.getState(), &cf_move.getState());
    }
}

TEST(Field, Interface)
{
    using FaceField = Block::Field<FaceData<double, AlignedBlockAllocator, 3>>;
    using IRange = typename FaceField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex faces(16);
    IRange face_domain(faces);
    FaceField ff(face_domain);
    EXPECT_EQ(ff.size(), face_domain.size());
    EXPECT_TRUE(ff.isScalar());
    EXPECT_NE(&ff.getState(), nullptr);
}

TEST(Field, Iterator)
{
    using CellField = Block::Field<CellData<float, AlignedBlockAllocator, 4>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(8);
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
    MIndex nodes(8);
    IRange node_domain(nodes);

    NodeField nf0(node_domain);
    FieldView nfv0(nf0);

    EXPECT_NE(&nf0, &nfv0); // separate instances

    EXPECT_EQ(nf0.getBlockPtr(), nfv0.getBlockPtr()); // same block data
    EXPECT_EQ(&nf0.getState(), &nfv0.getState());     // same state
    EXPECT_TRUE(nf0.isMemoryOwner());
    EXPECT_FALSE(nfv0.isMemoryOwner());

    // assignment with fields
    NodeField nf1(node_domain);
    std::fill(nf0.begin(), nf0.end(), 1);
    std::fill(nf1.begin(), nf1.end(), 2);
    int sum = 0;
    for (auto n : nf1) {
        sum += n;
    }
    EXPECT_EQ(sum, 2 * nf0.size());

    // 1. carrier = view
    EXPECT_NE(nf1.getBlockPtr(), nfv0.getBlockPtr());
    nf1 = nfv0; // deep copy
    EXPECT_NE(nf1.getBlockPtr(), nfv0.getBlockPtr());
    sum = 0;
    for (auto n : nf1) {
        sum += n;
    }
    EXPECT_EQ(sum, 1 * nf0.size());

    // 2. view = carrier
    EXPECT_NE(nf1.getBlockPtr(), nfv0.getBlockPtr());
    nfv0.setView(nf1); // shallow copy
    EXPECT_EQ(nf1.getBlockPtr(), nfv0.getBlockPtr());

    // 3. view = view
    FieldView nfv1(nf0);
    EXPECT_NE(nfv1.getBlockPtr(), nfv0.getBlockPtr());
    nfv1 = nfv0;
    EXPECT_EQ(nfv1.getBlockPtr(), nfv0.getBlockPtr());

    // 4. forced deep copies
    nfv0.setView(nf0);
    nfv1.setView(nf1);
    EXPECT_NE(nfv0.getBlockPtr(), nfv1.getBlockPtr());
    std::fill(nfv0.begin(), nfv0.end(), 1);
    std::fill(nfv1.begin(), nfv1.end(), 2);
    sum = 0;
    for (auto n : nfv1) {
        sum += n;
    }
    EXPECT_EQ(sum, 2 * nf0.size());

    // 4a. view = carrier
    nfv1.copyData(nf0);
    EXPECT_NE(nfv0.getBlockPtr(), nfv1.getBlockPtr());
    sum = 0;
    for (auto n : nfv1) {
        sum += n;
    }
    EXPECT_EQ(sum, 1 * nf0.size());

    std::fill(nfv1.begin(), nfv1.end(), 2);
    sum = 0;
    for (auto n : nfv1) {
        sum += n;
    }
    EXPECT_EQ(sum, 2 * nf0.size());

    // 4b. view = view
    nfv1.copyData(nfv0);
    EXPECT_NE(nfv0.getBlockPtr(), nfv1.getBlockPtr());
    sum = 0;
    for (auto n : nfv1) {
        sum += n;
    }
    EXPECT_EQ(sum, 1 * nf0.size());

    // 5. full copy view
    NodeField nf2 = nfv1.copy(); // call move constructor
    EXPECT_NE(nf2.getBlockPtr(), nfv1.getBlockPtr());
    sum = 0;
    for (auto n : nf2) {
        sum += n;
    }
    EXPECT_EQ(sum, 1 * nf0.size());

    std::fill(nfv0.begin(), nfv0.end(), 1);
    sum = 0;
    for (auto n : nfv0) {
        sum += n;
    }
    EXPECT_EQ(sum, 1 * nf0.size());
    std::fill(nfv1.begin(), nfv1.end(), 2);
    sum = 0;
    for (auto n : nfv1) {
        sum += n;
    }
    EXPECT_EQ(sum, 2 * nf0.size());
    EXPECT_NE(nfv0.getBlockPtr(), nfv1.getBlockPtr());
    nfv1.setView(nfv0);
    EXPECT_EQ(nfv0.getBlockPtr(), nfv1.getBlockPtr());
}

TEST(Field, Arithmetic)
{
    using CellField = Block::Field<CellData<float, AlignedBlockAllocator, 3>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using DataType = typename CellField::DataType;

    MIndex cells(8);
    IRange cell_domain(cells);
    CellField cf1(cell_domain);
    CellField cf2(cell_domain);
    const DataType s1 = 1;
    const DataType s2 = 2;
    std::fill(cf1.begin(), cf1.end(), s1);
    std::fill(cf2.begin(), cf2.end(), s2);

    // two field operands
    { // += | +
        auto cf(cf1);
        cf += cf2;
        for (const auto c : cf) {
            EXPECT_EQ(c, 3);
        }

        const auto cfp = cf1 + cf2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, 3);
        }
    }
    { // -= | -
        auto cf(cf1);
        cf -= cf2;
        for (const auto c : cf) {
            EXPECT_EQ(c, -1);
        }

        const auto cfp = cf1 - cf2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, -1);
        }
    }
    { // *= | *
        auto cf(cf1);
        cf *= cf2;
        for (const auto c : cf) {
            EXPECT_EQ(c, 2);
        }

        const auto cfp = cf1 * cf2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, 2);
        }
    }
    { // /= | /
        auto cf(cf1);
        cf /= cf2;
        for (const auto c : cf) {
            EXPECT_EQ(c, 0.5);
        }

        const auto cfp = cf1 / cf2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, 0.5);
        }
    }

    // one field operand
    { // += | +
        auto cf(cf1);
        cf += s2;
        for (const auto c : cf) {
            EXPECT_EQ(c, 3);
        }

        const auto cfp = cf1 + s2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, 3);
        }

        const auto cfpr = s2 + cf1;
        for (const auto c : cfpr) {
            EXPECT_EQ(c, 3);
        }
    }
    { // -= | -
        auto cf(cf1);
        cf -= s2;
        for (const auto c : cf) {
            EXPECT_EQ(c, -1);
        }

        const auto cfp = cf1 - s2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, -1);
        }

        const auto cfpr = s2 - cf1;
        for (const auto c : cfpr) {
            EXPECT_EQ(c, 1);
        }
    }
    { // *= | *
        auto cf(cf1);
        cf *= s2;
        for (const auto c : cf) {
            EXPECT_EQ(c, 2);
        }

        const auto cfp = cf1 * s2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, 2);
        }

        const auto cfpr = s2 * cf1;
        for (const auto c : cfpr) {
            EXPECT_EQ(c, 2);
        }
    }
    { // /= | /
        auto cf(cf1);
        cf /= s2;
        for (const auto c : cf) {
            EXPECT_EQ(c, 0.5);
        }

        const auto cfp = cf1 / s2;
        for (const auto c : cfp) {
            EXPECT_EQ(c, 0.5);
        }
    }

    // negation
    {
        const auto cf(-cf1);
        for (const auto c : cf) {
            EXPECT_EQ(c, -1);
        }
    }

    // reciprocal
    {
        auto cf(cf2);
        cf.reciprocal(2);
        for (const auto c : cf) {
            EXPECT_EQ(c, 1);
        }

        cf = cf2;
        cf.reciprocal();
        for (const auto c : cf) {
            EXPECT_EQ(c, 0.5);
        }
    }
}

TEST(FieldContainer, Construction)
{
    using NodeField = Block::Field<NodeData<char, AlignedBlockAllocator, 5>>;
    using IRange = typename NodeField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    using FV = Block::FieldView<NodeField>;
    using FC = Block::FieldContainer<NodeField>;
    MIndex nodes(8);
    IRange node_domain(nodes);

    { // default
        FC fc;
        EXPECT_EQ(fc.size(), 0);
    }

    { // construct new (owns memory)
        FC fc(2, node_domain);

        for (size_t i = 0; i < fc.size(); ++i) {
            NodeField &nf = fc[i];
            EXPECT_TRUE(nf.isMemoryOwner());
            std::fill(nf.begin(), nf.end(), 1);
        }
        int sum = 0;
        for (size_t i = 0; i < fc.size(); ++i) {
            for (const auto v : fc[i]) {
                sum += v;
            }
        }
        EXPECT_EQ(sum, 2 * node_domain.size());
    }

    { // construct from list of pointers (including FieldView types)
        FC fc(3, node_domain); // owns all
        FV fv0(fc[0]);         // view into first field

        std::vector<typename FC::FieldType *> ptr_list;
        ptr_list.push_back(&fv0);   // view
        ptr_list.push_back(&fc[1]); // owner
        ptr_list.push_back(&fc[2]); // owner
        FC fc1(ptr_list);

        EXPECT_EQ(fc1[0].getBlockPtr(), fc[0].getBlockPtr()); // view the same
        EXPECT_NE(fc1[1].getBlockPtr(), fc[1].getBlockPtr()); // copy
        EXPECT_NE(fc1[2].getBlockPtr(), fc[2].getBlockPtr()); // copy

        EXPECT_EQ(&fc1[0].getState(), &fc[0].getState()); // view the same
        EXPECT_NE(&fc1[1].getState(), &fc[1].getState()); // copy
        EXPECT_NE(&fc1[2].getState(), &fc[2].getState()); // copy
    }

    {                          // copy construction
        FC fc(4, node_domain); // owns all
        FV fv0(fc[0]);         // view
        FV fv1(fc[2]);         // view

        {                   // copy homogeneous
            FC fc_copy(fc); // deep copies
            for (size_t i = 0; i < fc_copy.size(); ++i) {
                EXPECT_TRUE(fc[i].isMemoryOwner());
                EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                EXPECT_NE(&fc_copy[i].getState(), &fc[i].getState());
                EXPECT_NE(&fc_copy[i], &fc[i]);
            }
        }
        { // copy mixed own/view
            std::vector<typename FC::FieldType *> ptr_list;
            ptr_list.push_back(&fv0);   // view
            ptr_list.push_back(&fc[1]); // owner
            ptr_list.push_back(&fv1);   // view
            ptr_list.push_back(&fc[3]); // owner
            FC fc1(ptr_list);
            FC fc_copy(fc1); // mixed
            for (size_t i = 0; i < fc_copy.size(); ++i) {
                if (i % 2 == 0) {
                    EXPECT_TRUE(fc[i].isMemoryOwner());
                    EXPECT_FALSE(fc_copy[i].isMemoryOwner());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                    EXPECT_EQ(&fc_copy[i].getState(), &fc1[i].getState());
                    EXPECT_EQ(&fc_copy[i].getState(), &fc[i].getState());
                } else {
                    EXPECT_TRUE(fc[i].isMemoryOwner());
                    EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                    EXPECT_NE(&fc_copy[i].getState(), &fc1[i].getState());
                    EXPECT_NE(&fc_copy[i].getState(), &fc[i].getState());
                }
                EXPECT_NE(&fc_copy[i], &fc[i]);
                EXPECT_NE(&fc_copy[i], &fc1[i]);
            }
        }
    }

    {                          // copy assignment
        FC fc(4, node_domain); // owns all
        FV fv0(fc[0]);         // view
        FV fv1(fc[2]);         // view

        {                 // copy homogeneous
            FC fc_copy;   // empty
            fc_copy = fc; // deep copies
            for (size_t i = 0; i < fc_copy.size(); ++i) {
                EXPECT_TRUE(fc[i].isMemoryOwner());
                EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                EXPECT_NE(&fc_copy[i].getState(), &fc[i].getState());
                EXPECT_NE(&fc_copy[i], &fc[i]);
            }
        }
        { // copy mixed own/view
            std::vector<typename FC::FieldType *> ptr_list;
            ptr_list.push_back(&fv0);   // view
            ptr_list.push_back(&fc[1]); // owner
            ptr_list.push_back(&fv1);   // view
            ptr_list.push_back(&fc[3]); // owner
            FC fc1(ptr_list);
            FC fc_copy;    // empty
            fc_copy = fc1; // mixed
            for (size_t i = 0; i < fc_copy.size(); ++i) {
                if (i % 2 == 0) {
                    EXPECT_TRUE(fc[i].isMemoryOwner());
                    EXPECT_FALSE(fc_copy[i].isMemoryOwner());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                    EXPECT_EQ(&fc_copy[i].getState(), &fc1[i].getState());
                    EXPECT_EQ(&fc_copy[i].getState(), &fc[i].getState());
                } else {
                    EXPECT_TRUE(fc[i].isMemoryOwner());
                    EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                    EXPECT_NE(&fc_copy[i].getState(), &fc1[i].getState());
                    EXPECT_NE(&fc_copy[i].getState(), &fc[i].getState());
                }
                EXPECT_NE(&fc_copy[i], &fc[i]);
                EXPECT_NE(&fc_copy[i], &fc1[i]);
            }
        }
        { // copy different sized containers
            std::vector<typename FC::FieldType *> ptr_list;
            ptr_list.push_back(&fc[1]); // owner
            ptr_list.push_back(&fv1);   // view
            ptr_list.push_back(&fc[3]); // owner
            FC fc1(ptr_list);           // 3 components
            EXPECT_EQ(fc1.size(), 3);
            fc1 = fc; // assign 4 components to a 3 component container (warns
                      // in debug build)
            EXPECT_EQ(fc1.size(), 4);
            for (size_t i = 0; i < ptr_list.size(); ++i) {
                EXPECT_NE(fc1[i].getBlockPtr(), ptr_list[i]);
                EXPECT_NE(&fc1[i].getState(), &fc[i].getState());
            }
        }
    }

    {                          // move construction
        FC fc(5, node_domain); // owns all
        const auto data = fc.getContainer();
        FC fc_move(std::move(fc));
        EXPECT_EQ(fc.size(), 0);
        EXPECT_EQ(fc_move.size(), data.size());
        for (size_t i = 0; i < fc_move.size(); ++i) {
            EXPECT_TRUE(fc_move[i].isMemoryOwner());
            EXPECT_EQ(fc_move[i].getBlockPtr(), data[i]->getBlockPtr());
            EXPECT_EQ(&fc_move[i].getState(), &data[i]->getState());
        }
    }

    {                          // move assignment
        FC fc(5, node_domain); // owns all
        const auto data = fc.getContainer();
        FC fc_move;
        fc_move = std::move(fc);
        EXPECT_EQ(fc.size(), 0);
        EXPECT_EQ(fc_move.size(), data.size());
        for (size_t i = 0; i < fc_move.size(); ++i) {
            EXPECT_TRUE(fc_move[i].isMemoryOwner());
            EXPECT_EQ(fc_move[i].getBlockPtr(), data[i]->getBlockPtr());
            EXPECT_EQ(&fc_move[i].getState(), &data[i]->getState());
        }
    }
}

TEST(FieldContainer, Iterator)
{
    using CellField = Block::Field<CellData<double, AlignedBlockAllocator, 5>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    using FC = Block::FieldContainer<CellField>;
    MIndex cells(8);
    IRange cell_domain(cells);
    FC field_container(8, cell_domain);

    // initialize
    double val = 0.0;
    for (auto block : field_container) {
        std::fill(block->begin(), block->end(), val);
        val += 1.0;
    }

    double sum = 0.0;
    double ref = 0.0;
    int count = 0;
    for (const auto block : field_container) {
        const auto &b = *block;
        ref += count * b.getBlockSize();
        ++count;
        for (const auto v : b) {
            sum += v;
        }
    }
    EXPECT_EQ(sum, ref);
}

TEST(FieldContainer, Interface)
{
    using NodeField = Block::Field<NodeData<size_t, AlignedBlockAllocator, 5>>;
    using IRange = typename NodeField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    using FC = Block::FieldContainer<NodeField>;
    MIndex nodes(16);
    IRange node_domain(nodes);

    {
        FC fc1;
        fc1.pushBack(nullptr);
        EXPECT_EQ(fc1.size(), 1);
        EXPECT_THROW(
            {
                try {
                    fc1[0];
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("FieldContainer: Component 0 was not assigned "
                                 "(nullptr)",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    fc1[0] = NodeField(node_domain);
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("FieldContainer: Component 0 was not assigned "
                                 "(nullptr)",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    fc1[Dir::X];
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("FieldContainer: Component 0 was not assigned "
                                 "(nullptr)",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
        EXPECT_THROW(
            {
                try {
                    fc1[Dir::X] = NodeField(node_domain);
                } catch (const std::runtime_error &e) {
                    EXPECT_STREQ("FieldContainer: Component 0 was not assigned "
                                 "(nullptr)",
                                 e.what());
                    throw;
                }
            },
            std::runtime_error);
    }

    {
        MIndex nodes(8);
        IRange node_domain(nodes);
        FC fc1(3, node_domain);
        FC fc2(3, node_domain);
        size_t k = 0;
        for (auto f : fc1) {
            std::fill(f->begin(), f->end(), k++);
        }
        fc2.copyData(fc1);
        k = 0;
        for (const auto f : fc2) {
            for (const auto n : *f) {
                EXPECT_EQ(n, k);
            }
            ++k;
        }
    }

    {
        MIndex nodes(8);
        IRange node_domain(nodes);
        FC fc;
        EXPECT_EQ(fc.size(), 0);
        for (size_t i = 1; i < 4; ++i) {
            fc.pushBack(new NodeField(node_domain));
            EXPECT_EQ(fc.size(), i);
        }
        EXPECT_EQ(fc.size(), 3);
        fc.clear();
        EXPECT_EQ(fc.size(), 0);
    }
}

TEST(FieldContainer, Arithmetic)
{
    using CellField = Block::Field<CellData<double, AlignedBlockAllocator, 2>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using DataType = typename CellField::DataType;

    using FC = Block::FieldContainer<CellField>;
    MIndex cells(8);
    IRange cell_domain(cells);
    FC cf1(3, cell_domain);
    FC cf2(3, cell_domain);
    const DataType s1 = 1;
    const DataType s2 = 2;
    for (auto f : cf1) {
        std::fill(f->begin(), f->end(), s1);
    }
    for (auto f : cf2) {
        std::fill(f->begin(), f->end(), s2);
    }

    // two field operands
    { // += | +
        auto cf(cf1);
        cf += cf2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 3);
            }
        }

        const auto cfp = cf1 + cf2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 3);
            }
        }
    }
    { // -= | -
        auto cf(cf1);
        cf -= cf2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, -1);
            }
        }

        const auto cfp = cf1 - cf2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, -1);
            }
        }
    }
    { // *= | *
        auto cf(cf1);
        cf *= cf2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 2);
            }
        }

        const auto cfp = cf1 * cf2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 2);
            }
        }
    }
    { // /= | /
        auto cf(cf1);
        cf /= cf2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 0.5);
            }
        }

        const auto cfp = cf1 / cf2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 0.5);
            }
        }
    }

    // one field operand
    { // += | +
        auto cf(cf1);
        cf += s2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 3);
            }
        }

        const auto cfp = cf1 + s2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 3);
            }
        }

        const auto cfpr = s2 + cf1;
        for (const auto f : cfpr) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 3);
            }
        }
    }
    { // -= | -
        auto cf(cf1);
        cf -= s2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, -1);
            }
        }

        const auto cfp = cf1 - s2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, -1);
            }
        }

        const auto cfpr = s2 - cf1;
        for (const auto f : cfpr) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 1);
            }
        }
    }
    { // *= | *
        auto cf(cf1);
        cf *= s2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 2);
            }
        }

        const auto cfp = cf1 * s2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 2);
            }
        }

        const auto cfpr = s2 * cf1;
        for (const auto f : cfpr) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 2);
            }
        }
    }
    { // /= | /
        auto cf(cf1);
        cf /= s2;
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 0.5);
            }
        }

        const auto cfp = cf1 / s2;
        for (const auto f : cfp) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 0.5);
            }
        }
    }

    // negation
    {
        const auto cf(-cf1);
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, -1);
            }
        }
    }

    // reciprocal
    {
        auto cf(cf2);
        cf.reciprocal(2);
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 1);
            }
        }

        cf = cf2;
        cf.reciprocal();
        for (const auto f : cf) {
            for (const auto c : *f) {
                EXPECT_EQ(c, 0.5);
            }
        }
    }
}

TEST(TensorField, Construction)
{
    // Rank 2 tensor, cell-centered data, 3-dimensional
    using TensorField =
        Block::TensorField<double, 2, Cubism::EntityType::Cell, 3>;
    using IRange = typename TensorField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    TensorField tf(cell_domain);

    // If you like that better; must be castable to size_t and 0 for first index
    enum MyIndex { XX = 0, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };
    EXPECT_EQ(tf[MyIndex::XY].getBlockSize(), cell_domain.size());
    EXPECT_EQ(tf.size(), std::pow(IRange::Dim, TensorField::Rank));
    const auto *pstate = &tf.getState();
    for (auto c : tf) {
        EXPECT_EQ(&c->getState(), pstate);
    }

    { // low-level
        TensorField tf_copy(tf,
                            TensorField::FieldType::BaseType::MemoryOwner::Yes);
        EXPECT_EQ(tf_copy.size(), tf.size());
        const auto *pstate = &tf_copy.getState();
        for (auto c : tf_copy) {
            EXPECT_EQ(&c->getState(), pstate);
        }
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tf[i].getBlockPtr(), tf_copy[i].getBlockPtr());
            EXPECT_NE(&tf[i].getState(), &tf_copy[i].getState());
        }
    }

    { // copy construction
        TensorField tf1(tf);
        EXPECT_EQ(tf1.size(), tf.size());
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tf[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tf[i].getState(), &tf1[i].getState());
        }
    }

    { // copy assignment
        TensorField tf1(tf);
        TensorField tf2(cell_domain);
        tf1 = tf2;
        EXPECT_EQ(tf1.size(), tf.size());
        const auto *pstate = &tf1.getState();
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tf[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(tf1[i].getBlockPtr(), tf2[i].getBlockPtr());
            EXPECT_NE(tf[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tf[i].getState(), &tf1[i].getState());
            EXPECT_EQ(pstate, &tf1[i].getState());
        }
    }

    { // move construction
        using FieldView = Block::FieldView<TensorField>;
        TensorField tf_copy(tf);
        FieldView tv(tf);
        TensorField tf1(std::move(tf));
        EXPECT_EQ(tf.size(), 0);
        EXPECT_EQ(tf1.size(), tv.size());
        const auto *pstate = &tf1.getState();
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tv[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tv[i].getState(), &tf1[i].getState());
            EXPECT_EQ(pstate, &tf1[i].getState());
        }
        tf = tf_copy;
    }

    { // move assignment
        using FieldView = Block::FieldView<TensorField>;
        TensorField tf_copy(tf);
        FieldView tv(tf);
        TensorField tf1(cell_domain);
        tf1 = std::move(tf);
        EXPECT_EQ(tf.size(), 0);
        EXPECT_EQ(tf1.size(), tv.size());
        const auto *pstate = &tf1.getState();
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tv[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tv[i].getState(), &tf1[i].getState());
            EXPECT_EQ(pstate, &tf1[i].getState());
        }
        tf = tf_copy;
    }
}

TEST(FaceContainer, Construction)
{
    using FaceField = Block::FaceContainer<
        Block::Field<FaceData<double, AlignedBlockAllocator, 3>>>;
    using IRange = typename FaceField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    FaceField ff(cell_domain);
    EXPECT_EQ(ff.size(), CUBISM_DIMENSION);
    const auto *pstate = &ff.getState();
    for (auto c : ff) {
        EXPECT_EQ(&c->getState(), pstate);
    }

    { // low-level
        FaceField ff_copy(ff, FaceField::FieldType::BaseType::MemoryOwner::Yes);
        EXPECT_EQ(ff_copy.size(), ff.size());
        const auto *pstate = &ff_copy.getState();
        for (auto c : ff_copy) {
            EXPECT_EQ(&c->getState(), pstate);
        }
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(ff[i].getBlockPtr(), ff_copy[i].getBlockPtr());
            EXPECT_NE(&ff[i].getState(), &ff_copy[i].getState());
        }
    }

    { // copy construction
        FaceField ff1(ff);
        EXPECT_EQ(ff1.size(), ff.size());
        const auto *pstate = &ff1.getState();
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(ff[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&ff[i].getState(), &ff1[i].getState());
            EXPECT_EQ(pstate, &ff1[i].getState());
        }
    }

    { // copy assignment
        FaceField ff1(ff);
        FaceField ff2(cell_domain);
        ff1 = ff2;
        EXPECT_EQ(ff1.size(), ff.size());
        EXPECT_EQ(ff1.size(), ff2.size());
        const auto *pstate = &ff1.getState();
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(ff[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(ff1[i].getBlockPtr(), ff2[i].getBlockPtr());
            EXPECT_NE(&ff[i].getState(), &ff1[i].getState());
            EXPECT_NE(&ff1[i].getState(), &ff2[i].getState());
            EXPECT_EQ(pstate, &ff1[i].getState());
        }
    }

    { // move construction
        using FieldView = Block::FieldView<FaceField>;
        FaceField ff_copy(ff);
        FieldView fv(ff);
        FaceField ff1(std::move(ff));
        EXPECT_EQ(ff.size(), 0);
        EXPECT_EQ(ff1.size(), fv.size());
        const auto *pstate = &ff1.getState();
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(fv[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&fv[i].getState(), &ff1[i].getState());
            EXPECT_EQ(pstate, &ff1[i].getState());
        }

        ff = ff_copy;
    }

    { // move assignment
        using FieldView = Block::FieldView<FaceField>;
        FaceField ff_copy(ff);
        FieldView fv(ff);
        FaceField ff1(cell_domain);
        ff1 = std::move(ff);
        EXPECT_EQ(ff.size(), 0);
        EXPECT_EQ(ff1.size(), fv.size());
        const auto *pstate = &ff1.getState();
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(fv[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&fv[i].getState(), &ff1[i].getState());
            EXPECT_EQ(pstate, &ff1[i].getState());
        }

        ff = ff_copy;
    }
}

TEST(FaceContainer, TensorConstruction)
{
    using TensorField =
        Block::TensorField<double, 2, Cubism::EntityType::Face, 3>;
    using FaceField = Block::FaceContainer<TensorField>;
    using IRange = typename FaceField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    FaceField ff(cell_domain);
    EXPECT_EQ(ff.size(), CUBISM_DIMENSION);
    const auto *pstate = &ff.getState();
    for (auto f : ff) { // individual face directions
        EXPECT_EQ(&f->getState(), pstate);
        for (auto c : *f) { // tensor components for each f
            EXPECT_EQ(&c->getState(), pstate);
        }
    }
}

TEST(FieldView, Construction)
{
    using CellField = Block::Field<CellData<double, AlignedBlockAllocator, 3>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using FieldView = Block::FieldView<CellField>;

    MIndex cells(16);
    IRange cell_domain(cells);
    CellField cf(cell_domain);
    FieldView cf_view(cf);
    EXPECT_TRUE(cf.isMemoryOwner());
    EXPECT_FALSE(cf_view.isMemoryOwner());
    EXPECT_EQ(cf.getBlockPtr(), cf_view.getBlockPtr());
    EXPECT_EQ(&cf.getState(), &cf_view.getState());
}

TEST(FieldView, Copy)
{
    using CellField = Block::Field<CellData<double, AlignedBlockAllocator, 3>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using FieldView = Block::FieldView<CellField>;

    MIndex cells(16);
    IRange cell_domain(cells);
    CellField cf(cell_domain);
    FieldView cf_view(cf);
    EXPECT_TRUE(cf.isMemoryOwner());
    EXPECT_FALSE(cf_view.isMemoryOwner());

    { // copy constructor
        typename FieldView::BaseType cf_copy = cf_view.copy();
        EXPECT_TRUE(cf_copy.isMemoryOwner());
        EXPECT_NE(cf.getBlockPtr(), cf_copy.getBlockPtr());
        EXPECT_NE(&cf.getState(), &cf_copy.getState());
    }

    { // copy assignment (a move in this case)
        CellField cf_copy(cell_domain);

        cf_copy = cf_view.copy();
        EXPECT_TRUE(cf_copy.isMemoryOwner());
        EXPECT_NE(cf.getBlockPtr(), cf_copy.getBlockPtr());
        EXPECT_NE(&cf.getState(), &cf_copy.getState());
    }
}

} // namespace
