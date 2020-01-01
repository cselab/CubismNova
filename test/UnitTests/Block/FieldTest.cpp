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

TEST(FieldContainer, Construction)
{
    using NodeField = Block::Field<NodeData<char, AlignedBlockAllocator, 5>>;
    using FieldState = typename NodeField::FieldStateType;
    using IRange = typename NodeField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    using FV = Block::FieldView<NodeField>;
    using FC = Block::FieldContainer<NodeField>;
    MIndex nodes(16);
    IRange node_domain(nodes);

    { // default
        FC fc;
        EXPECT_EQ(fc.size(), 0);

        FC fc1(1);
        EXPECT_EQ(fc1.size(), 1);
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
        FV fv0(fc[0]);

        std::vector<typename FC::FieldType *> ptr_list;
        ptr_list.push_back(&fv0);   // view
        ptr_list.push_back(&fc[1]); // owner
        ptr_list.push_back(&fc[2]); // owner
        FC fc1(ptr_list);

        EXPECT_EQ(fc1[0].getBlockPtr(), fc[0].getBlockPtr()); // the same
        EXPECT_NE(fc1[1].getBlockPtr(), fc[1].getBlockPtr()); // copy
        EXPECT_NE(fc1[2].getBlockPtr(), fc[2].getBlockPtr()); // copy

        // view of a field container
        using FVC = Block::FieldView<FC>;
        FVC fvc(fc1);
        EXPECT_FALSE(fvc[0].isMemoryOwner());
        EXPECT_FALSE(fvc[1].isMemoryOwner());
        EXPECT_FALSE(fvc[2].isMemoryOwner());
        EXPECT_EQ(fc1[0].getBlockPtr(), fvc[0].getBlockPtr());
        EXPECT_EQ(fc1[1].getBlockPtr(), fvc[1].getBlockPtr());
        EXPECT_EQ(fc1[2].getBlockPtr(), fvc[2].getBlockPtr());
        EXPECT_NE(&fc1[0], &fvc[0]);
        EXPECT_NE(&fc1[1], &fvc[1]);
        EXPECT_NE(&fc1[2], &fvc[2]);
    }

    // low-level constructors
    {
        FC fc(2, node_domain);
        std::vector<IRange> rl;
        std::vector<typename NodeField::DataType *> pl;
        std::vector<size_t> bl;
        std::vector<FieldState *> sl;
        for (size_t i = 0; i < fc.size(); ++i) {
            rl.push_back(fc[i].getIndexRange());
            pl.push_back(fc[i].getData());
            bl.push_back(fc[i].getBlockBytes());
            sl.push_back(&fc[i].getState());
        }
        FC fc1(rl, pl, bl, sl); // never owns data
        for (size_t i = 0; i < fc.size(); ++i) {
            EXPECT_EQ(fc1[i].getRank(), 0);
            EXPECT_EQ(fc1[i].getComp(), i);
            EXPECT_FALSE(fc1[i].isMemoryOwner());
            EXPECT_EQ(fc1[i].getBlockPtr(), fc[i].getBlockPtr());
            EXPECT_NE(&fc1[i], &fc[i]);
        }
    }

    {                          // copy construction
        FC fc(4, node_domain); // owns all
        FV fv0(fc[0]);         // view
        FV fv1(fc[2]);         // view

        {                   // copy homogeneous
            FC fc_copy(fc); // deep copies
            for (size_t i = 0; i < fc_copy.size(); ++i) {
                EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
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
                    EXPECT_FALSE(fc_copy[i].isMemoryOwner());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                } else {
                    EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
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
                EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
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
                    EXPECT_FALSE(fc_copy[i].isMemoryOwner());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_EQ(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                } else {
                    EXPECT_TRUE(fc_copy[i].isMemoryOwner());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc1[i].getBlockPtr());
                    EXPECT_NE(fc_copy[i].getBlockPtr(), fc[i].getBlockPtr());
                }
                EXPECT_NE(&fc_copy[i], &fc[i]);
                EXPECT_NE(&fc_copy[i], &fc1[i]);
            }
        }
    }

    {                                     // move construction
        FC fc(5, node_domain);            // owns all
        using FVC = Block::FieldView<FC>; // view of a field container
        FVC fvc(fc);
        FC fc_move(std::move(fc));
        EXPECT_EQ(fc.size(), 0);
        EXPECT_EQ(fc_move.size(), fvc.size());
        for (size_t i = 0; i < fc_move.size(); ++i) {
            EXPECT_TRUE(fc_move[i].isMemoryOwner());
            EXPECT_FALSE(fvc[i].isMemoryOwner());
            EXPECT_EQ(fc_move[i].getBlockPtr(), fvc[i].getBlockPtr());
            EXPECT_NE(&fc_move[i], &fvc[i]);
        }
    }

    {                                     // move assignment
        FC fc(5, node_domain);            // owns all
        using FVC = Block::FieldView<FC>; // view of a field container
        FVC fvc(fc);
        FC fc_move;
        fc_move = std::move(fc);
        EXPECT_EQ(fc.size(), 0);
        EXPECT_EQ(fc_move.size(), fvc.size());
        for (size_t i = 0; i < fc_move.size(); ++i) {
            EXPECT_TRUE(fc_move[i].isMemoryOwner());
            EXPECT_FALSE(fvc[i].isMemoryOwner());
            EXPECT_EQ(fc_move[i].getBlockPtr(), fvc[i].getBlockPtr());
            EXPECT_NE(&fc_move[i], &fvc[i]);
        }
    }
}

TEST(FieldContainer, Iterator)
{
    using CellField = Block::Field<CellData<double, AlignedBlockAllocator, 5>>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    using FC = Block::FieldContainer<CellField>;
    MIndex cells(16);
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
    for (auto block : field_container) {
        auto &b = *block;
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
        FC fc1(1);
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
        for (auto f : fc2) {
            for (const auto n : *f) {
                EXPECT_EQ(n, k);
            }
            ++k;
        }
    }
}

// TODO: [fabianw@mavt.ethz.ch; 2020-01-01]
TEST(FieldContainer, Arithmetic) {}

TEST(FaceFieldAll, Construction)
{
    // CUBISM_DIMENSION-ional FaceField
    using FaceField = Block::FaceFieldAll<double>;
    using IRange = typename FaceField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    FaceField ff(cell_domain);
    EXPECT_EQ(ff.size(), CUBISM_DIMENSION);
    size_t k = 0;
    for (auto f : ff) {
        EXPECT_EQ(f->getRank(), 0);
        EXPECT_EQ(f->getComp(), k++);
    }

    { // low-level
        std::vector<IRange> rl;
        std::vector<typename FaceField::DataType *> pl;
        std::vector<size_t> bl;
        std::vector<typename FaceField::FieldState *> sl;

        FaceField ff_copy(ff, FaceField::FieldType::BaseType::MemoryOwner::Yes);
        EXPECT_EQ(ff_copy.size(), ff.size());
        for (size_t i = 0; i < ff.size(); ++i) {
            rl.push_back(ff[i].getIndexRange());
            pl.push_back(ff[i].getData());
            bl.push_back(ff[i].getBlockBytes());
            sl.push_back(&ff[i].getState());
            EXPECT_NE(ff[i].getBlockPtr(), ff_copy[i].getBlockPtr());
            EXPECT_NE(&ff[i].getState(), &ff_copy[i].getState());
        }

        FaceField ff_view(rl, pl, bl, sl);
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_EQ(ff[i].getBlockPtr(), ff_view[i].getBlockPtr());
            EXPECT_EQ(&ff[i].getState(), &ff_view[i].getState());
        }
    }

    { // copy construction
        FaceField ff1(ff);
        EXPECT_EQ(ff1.size(), ff.size());
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(ff[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&ff[i].getState(), &ff1[i].getState());
        }
    }

    { // copy assignment
        FaceField ff1;
        EXPECT_EQ(ff1.size(), 0);
        ff1 = ff;
        EXPECT_EQ(ff1.size(), ff.size());
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(ff[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&ff[i].getState(), &ff1[i].getState());
        }
    }

    { // move construction
        using FieldView = Block::FieldView<FaceField>;
        FaceField ff_copy(ff);
        FieldView fv(ff);
        FaceField ff1(std::move(ff));
        EXPECT_EQ(ff.size(), 0);
        EXPECT_EQ(ff1.size(), fv.size());
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(fv[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&fv[i].getState(), &ff1[i].getState());
        }

        ff = ff_copy;
    }

    { // move assignment
        using FieldView = Block::FieldView<FaceField>;
        FaceField ff_copy(ff);
        FieldView fv(ff);
        FaceField ff1;
        ff1 = std::move(ff);
        EXPECT_EQ(ff.size(), 0);
        EXPECT_EQ(ff1.size(), fv.size());
        for (size_t i = 0; i < ff.size(); ++i) {
            EXPECT_NE(fv[i].getBlockPtr(), ff1[i].getBlockPtr());
            EXPECT_NE(&fv[i].getState(), &ff1[i].getState());
        }

        ff = ff_copy;
    }
}

TEST(TensorField, Construction)
{
    using CellField = Block::Field<CellData<double, AlignedBlockAllocator, 3>>;
    using TensorField = Block::TensorField<CellField, 2>; // Rank 2 tensor
    using IRange = typename TensorField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);
    IRange cell_domain(cells);
    TensorField tf(cell_domain);

    // If you like that better; must be castable to size_t and 0 for first index
    enum MyIndex { XX = 0, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };
    EXPECT_EQ(tf[MyIndex::XY].getBlockSize(), cell_domain.size());

    EXPECT_EQ(tf.size(), std::pow(IRange::Dim, TensorField::Rank));
    size_t k = 0;
    for (auto c : tf) {
        EXPECT_EQ(c->getRank(), TensorField::Rank);
        EXPECT_EQ(c->getComp(), k++);
    }

    { // low-level
        std::vector<IRange> rl;
        std::vector<typename CellField::DataType *> pl;
        std::vector<size_t> bl;
        std::vector<typename TensorField::FieldState *> sl;

        TensorField tf_copy(tf,
                            TensorField::FieldType::BaseType::MemoryOwner::Yes);
        EXPECT_EQ(tf_copy.size(), tf.size());
        for (size_t i = 0; i < tf.size(); ++i) {
            rl.push_back(tf[i].getIndexRange());
            pl.push_back(tf[i].getData());
            bl.push_back(tf[i].getBlockBytes());
            sl.push_back(&tf[i].getState());
            EXPECT_NE(tf[i].getBlockPtr(), tf_copy[i].getBlockPtr());
            EXPECT_NE(&tf[i].getState(), &tf_copy[i].getState());
        }

        TensorField tf_view(rl, pl, bl, sl);
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_EQ(tf[i].getBlockPtr(), tf_view[i].getBlockPtr());
            EXPECT_EQ(&tf[i].getState(), &tf_view[i].getState());
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
        TensorField tf1;
        EXPECT_EQ(tf1.size(), 0);
        tf1 = tf;
        EXPECT_EQ(tf1.size(), tf.size());
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tf[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tf[i].getState(), &tf1[i].getState());
        }
    }

    { // move construction
        using FieldView = Block::FieldView<TensorField>;
        TensorField tf_copy(tf);
        FieldView tv(tf);
        TensorField tf1(std::move(tf));
        EXPECT_EQ(tf.size(), 0);
        EXPECT_EQ(tf1.size(), tv.size());
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tv[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tv[i].getState(), &tf1[i].getState());
        }

        tf = tf_copy;
    }

    { // move assignment
        using FieldView = Block::FieldView<TensorField>;
        TensorField tf_copy(tf);
        FieldView tv(tf);
        TensorField tf1;
        tf1 = std::move(tf);
        EXPECT_EQ(tf.size(), 0);
        EXPECT_EQ(tf1.size(), tv.size());
        for (size_t i = 0; i < tf.size(); ++i) {
            EXPECT_NE(tv[i].getBlockPtr(), tf1[i].getBlockPtr());
            EXPECT_NE(&tv[i].getState(), &tf1[i].getState());
        }

        tf = tf_copy;
    }
}

} // namespace
