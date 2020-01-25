// File       : FieldAOSTest.cpp
// Created    : Tue Jan 21 2020 06:25:58 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Convert SoA fields to AoS for I/O
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "IO/FieldAOS.h"
#include "Block/Field.h"
#include "gtest/gtest.h"
#include <algorithm>

namespace
{
using namespace Cubism;

TEST(IO, FieldAOS)
{
    using CellField = Block::CellField<float>;
    using DataType = typename CellField::DataType;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using VecCellField = Block::VectorField<float, Cubism::EntityType::Cell>;

    MIndex cells(8);
    IRange cell_domain(cells);

    {
        // scalar
        CellField cf(cell_domain);
        std::fill(cf.begin(), cf.end(), 1.234);
        DataType *buf = new DataType[cell_domain.size()];
        IO::Field2AOS(cf, cell_domain, buf);
        for (size_t i = 0; i < cf.size(); ++i) {
            EXPECT_EQ(buf[i], cf[i]);
        }
        delete[] buf;

        // vector
        VecCellField vf(cell_domain);
        int k = 1;
        for (auto c : vf) {
            std::fill(c->begin(), c->end(), k * 1.234);
            k += 1;
        }
        buf = new DataType[VecCellField::NComponents * cell_domain.size()];
        IO::Field2AOS(vf, cell_domain, buf);
        for (size_t i = 0; i < cf.size(); ++i) {
            for (size_t c = 0; c < VecCellField::NComponents; ++c) {
                EXPECT_EQ(buf[c + i * VecCellField::NComponents], vf[c][i]);
            }
        }
        delete[] buf;
    }
    { // scalar subspace
        CellField cf(cell_domain);
        std::fill(cf.begin(), cf.end(), 0);
        IRange subrange(2, 5); // 5 is exclusive
        for (auto &p : cf.getIndexRange()) {
            if (subrange.isGlobalIndex(p)) {
                cf[p] = 1.234;
            }
        }
        DataType *buf = new DataType[subrange.size()];
        IO::Field2AOS(cf, subrange, buf);
        for (auto &p : cf.getIndexRange()) {
            if (subrange.isGlobalIndex(p)) {
                EXPECT_EQ(buf[subrange.getFlatIndexFromGlobal(p)], cf[p]);
            } else {
                EXPECT_EQ(0, cf[p]);
            }
        }
        delete[] buf;
    }
    { // vector subspace
        VecCellField vf(cell_domain);
        int k = 1;
        for (auto c : vf) {
            std::fill(c->begin(), c->end(), k * 0.5);
            k += 1;
        }

        k = 1;
        IRange subrange(2, 5); // 5 is exclusive
        for (auto c : vf) {
            auto &bf = *c;
            for (auto &p : bf.getIndexRange()) {
                if (subrange.isGlobalIndex(p)) {
                    bf[p] = k;
                }
            }
            k += 1;
        }
        DataType *buf =
            new DataType[VecCellField::NComponents * subrange.size()];
        IO::Field2AOS(vf, subrange, buf);
        k = 1;
        for (auto c : vf) {
            auto &bf = *c;
            for (auto &p : bf.getIndexRange()) {
                if (subrange.isGlobalIndex(p)) {
                    EXPECT_EQ(
                        buf[(k - 1) + VecCellField::NComponents *
                                          subrange.getFlatIndexFromGlobal(p)],
                        k);
                } else {
                    EXPECT_EQ(k * 0.5, bf[p]);
                }
            }
            k += 1;
        }
        delete[] buf;
    }
}
} // namespace
