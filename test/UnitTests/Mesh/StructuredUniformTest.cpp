// File       : StructuredUniformTest.cpp
// Created    : Fri Jan 03 2020 03:19:05 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Structured uniform mesh test
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Core/Index.h"
#include "gtest/gtest.h"
#include <iostream>

template class Cubism::Mesh::StructuredUniform<double, 3>;
template class Cubism::Mesh::StructuredUniform<double, 4>;
template class Cubism::Mesh::StructuredUniform<float, 2>;
template class Cubism::Mesh::StructuredUniform<float, 3>;

namespace
{
using namespace Cubism;

TEST(StructuredUniform, Construction)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;
    using Entity = typename Mesh::EntityType;
    using Range = typename Mesh::RangeType;

    std::cout << "sizeof(Mesh::StructuredUniform) = " << sizeof(Mesh)
              << std::endl;

    const PointType start(-1);
    const PointType end(1);
    const PointType extent = end - start;
    const MIndex cells{6, 7, 8};
    const PointType h = extent / PointType(cells);

    { // local and global origin at 0
        Mesh m(end, cells, MeshIntegrity::FullMesh);
        EXPECT_EQ(m.getExtent(), extent / 2);
        EXPECT_EQ(m.getOrigin(), MIndex(0));
        EXPECT_EQ(m.getGlobalOrigin(), MIndex(0));
        EXPECT_EQ(m.getMultiIndex(0, Entity::Cell), MIndex(0));
        for (const auto &c : m.getIterator(Entity::Cell)) {
            EXPECT_EQ(m.getCellSize(c), h / 2);
        }
    }

    { // local and global origin at start
        Mesh m(start, end, cells, MeshIntegrity::FullMesh);
        EXPECT_EQ(m.getExtent(), extent);
        EXPECT_EQ(m.getOrigin(), start);
        EXPECT_EQ(m.getGlobalOrigin(), start);
        EXPECT_EQ(m.getMultiIndex(0, Entity::Cell), MIndex(0));
        for (const auto &c : m.getIterator(Entity::Cell)) {
            EXPECT_EQ(m.getCellSize(c), h);
        }
    }

    { // local origin at 0; global origin at start
        const PointType gorigin = start;
        const Range phys_domain(MIndex(0), end);
        const IRange cell_domain(cells);
        Mesh m(gorigin, phys_domain, cell_domain, MeshIntegrity::FullMesh);
        EXPECT_EQ(m.getExtent(), extent / 2);
        EXPECT_EQ(m.getOrigin(), MIndex(0));
        EXPECT_EQ(m.getGlobalOrigin(), start);
        EXPECT_EQ(m.getMultiIndex(0, Entity::Cell), MIndex(0));
        for (const auto &c : m.getIterator(Entity::Cell)) {
            EXPECT_EQ(m.getCellSize(c), h / 2);
        }
    }

    { // low-level
        const PointType gorigin = start;
        const Range phys_domain(MIndex(0), end);
        Mesh m(
            gorigin,
            phys_domain,
            IRange(cells),                                 // cells
            IRange(cells + 1),                             // nodes
            std::vector<IRange>(Mesh::Dim, IRange(cells)), // faces base domain
            MeshIntegrity::FullMesh);
        EXPECT_EQ(m.getExtent(), extent / 2);
        EXPECT_EQ(m.getOrigin(), MIndex(0));
        EXPECT_EQ(m.getGlobalOrigin(), start);
        EXPECT_EQ(m.getMultiIndex(0, Entity::Cell), MIndex(0));
        for (const auto &c : m.getIterator(Entity::Cell)) {
            EXPECT_EQ(m.getCellSize(c), h / 2);
        }
    }

    {
        Mesh m0(end, cells, MeshIntegrity::FullMesh);
        Mesh m1(m0);
        for (const auto &c : m0.getIterator(Entity::Cell)) {
            EXPECT_EQ(m0.getCellSize(c), m1.getCellSize(c));
        }
    }

    {
        Mesh m0(end, cells, MeshIntegrity::FullMesh);
        Mesh m1(m0);
        Mesh m2(std::move(m1));
        for (const auto &c : m0.getIterator(Entity::Cell)) {
            EXPECT_EQ(m0.getCellSize(c), m2.getCellSize(c));
        }
    }
}

TEST(StructuredUniform, Iterator)
{
    using IRange = Core::IndexRange<4>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;
    using Entity = typename Mesh::EntityType;

    const PointType end(1);
    const MIndex cells(4);
    const PointType h = 1 / PointType(cells);
    Mesh m(end, cells, MeshIntegrity::FullMesh);

    { // cells
        for (const auto &c : m.getIterator(Entity::Cell)) {
            EXPECT_EQ(m.getCellSize(c), h);
        }
    }
    { // nodes
        for (const auto &n : m.getIterator(Entity::Node)) {
            EXPECT_EQ(m.getCellSize(n),
                      h); // uniform cell size applies to all entities
        }
    }
    { // faces (4-dimensional domain)
        for (size_t d = 0; d < 4; ++d) {
            for (const auto &f : m.getIterator(Entity::Face, d)) {
                EXPECT_EQ(m.getCellSize(f),
                          h); // uniform cell size applies to all entities
            }
        }
    }
    { // faces (4-dimensional domain, custom direction indicator)
        enum class Dir { X1 = 0, X2, X3, X4 };
        for (const auto &f :
             m.getIterator(Entity::Face, Dir::X1)) { // direction 1
            EXPECT_EQ(m.getCellSize(f),
                      h); // uniform cell size applies to all entities
        }
        for (const auto &f :
             m.getIterator(Entity::Face, Dir::X2)) { // direction 2
            EXPECT_EQ(m.getCellSize(f),
                      h); // uniform cell size applies to all entities
        }
        for (const auto &f :
             m.getIterator(Entity::Face, Dir::X3)) { // direction 3
            EXPECT_EQ(m.getCellSize(f),
                      h); // uniform cell size applies to all entities
        }
        for (const auto &f :
             m.getIterator(Entity::Face, Dir::X4)) { // direction 4
            EXPECT_EQ(m.getCellSize(f),
                      h); // uniform cell size applies to all entities
        }
    }
}

TEST(StructuredUniform, BasicInterface)
{
    using IRange = Core::IndexRange<2>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<float, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using Range = typename Mesh::RangeType;
    using RealType = typename Mesh::RealType;
    using PointType = typename Mesh::PointType;
    using Entity = typename Mesh::EntityType;

    const PointType origin(-2);
    const PointType start{-1.333, -0.6789};
    const PointType end(1);
    const MIndex cells{4, 2};
    const PointType domain = end - start;
    const PointType h = (end - start) / PointType(cells);
    Mesh m(origin, Range(start, end), IRange(cells), MeshIntegrity::FullMesh);

    { // iterator cell
        auto it = m.getIterator(Entity::Cell).begin();
        auto ite = m.getIterator(Entity::Cell).end();
        auto it0(it);
        const MIndex i0(0);
        const MIndex ix = MIndex::getUnitVector(0);
        const MIndex iy = MIndex::getUnitVector(1);
        EXPECT_EQ(it, it0);
        EXPECT_EQ(*it++, i0 + 0 * ix + 0 * iy);
        EXPECT_NE(it, it0);
        it0 = it;
        EXPECT_EQ(it, it0);
        EXPECT_EQ(*it++, i0 + 1 * ix + 0 * iy);
        EXPECT_NE(it, it0);
        EXPECT_EQ(*it++, i0 + 2 * ix + 0 * iy);
        EXPECT_EQ(*it++, i0 + 3 * ix + 0 * iy);
        EXPECT_EQ(*it++, i0 + 0 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 1 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 2 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 3 * ix + 1 * iy);
        EXPECT_EQ(it, ite);
    }
    { // iterator face-x
        auto it = m.getIterator(Entity::Face, Dir::X).begin();
        auto ite = m.getIterator(Entity::Face, Dir::X).end();
        auto it0(it);
        const MIndex i0(0);
        const MIndex ix = MIndex::getUnitVector(Dir::X);
        const MIndex iy = MIndex::getUnitVector(Dir::Y);
        EXPECT_EQ(it, it0);
        EXPECT_EQ(*it++, i0 + 0 * ix + 0 * iy);
        EXPECT_NE(it, it0);
        it0 = it;
        EXPECT_EQ(it, it0);
        EXPECT_EQ(*it++, i0 + 1 * ix + 0 * iy);
        EXPECT_NE(it, it0);
        EXPECT_EQ(*it++, i0 + 2 * ix + 0 * iy);
        EXPECT_EQ(*it++, i0 + 3 * ix + 0 * iy);
        EXPECT_EQ(*it++, i0 + 4 * ix + 0 * iy);
        EXPECT_EQ(*it++, i0 + 0 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 1 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 2 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 3 * ix + 1 * iy);
        EXPECT_EQ(*it++, i0 + 4 * ix + 1 * iy);
        EXPECT_EQ(it, ite);
    }
    { // high-level mesh iteration
        // cells
        for (const auto &c : m[Entity::Cell]) {
            const auto p = m.getGlobalIndex(c, Entity::Cell);
            EXPECT_EQ(p, c);
        }

        // nodes
        for (const auto &n : m[Entity::Node]) {
            const auto p = m.getGlobalIndex(n, Entity::Node);
            EXPECT_EQ(p, n);
        }

        // faces (variants)
        for (const auto &f : m[Entity::Face][Dir::Y]) {
            const auto p = m.getGlobalIndex(f, Entity::Face, Dir::Y);
            EXPECT_EQ(p, f);
        }
        for (size_t d = 0; d < Mesh::Dim; ++d) {
            for (const auto &f : m[Entity::Face][d]) {
                const auto p = m.getGlobalIndex(f, Entity::Face, d);
                EXPECT_EQ(p, f);
            }
        }
    }
    { // sizes
        EXPECT_EQ(m.size(Entity::Cell), 4 * 2);
        EXPECT_EQ(m.size(Entity::Node), (4 + 1) * (2 + 1));
        EXPECT_EQ(m.size(Entity::Face, Dir::X), (4 + 1) * 2);
        EXPECT_EQ(m.size(Entity::Face, Dir::Y), 4 * (2 + 1));

        MIndex p0{4, 2};
        EXPECT_EQ(m.getSize(Entity::Cell), p0);
        p0 = MIndex{4 + 1, 2 + 1};
        EXPECT_EQ(m.getSize(Entity::Node), p0);
        p0 = MIndex{4 + 1, 2};
        EXPECT_EQ(m.getSize(Entity::Face, Dir::X), p0);
        p0 = MIndex{4, 2 + 1};
        EXPECT_EQ(m.getSize(Entity::Face, Dir::Y), p0);
    }
    { // index ranges
        MIndex p0(0);
        EXPECT_EQ(m.getIndexRange(Entity::Cell).getBegin(), p0);
        EXPECT_EQ(m.getIndexRange(Entity::Cell).getEnd(), cells);
        EXPECT_EQ(m.getIndexRange(Entity::Node).getBegin(), p0);
        EXPECT_EQ(m.getIndexRange(Entity::Node).getEnd(), cells + 1);
        EXPECT_EQ(m.getIndexRange(Entity::Face, Dir::X).getBegin(), p0);
        EXPECT_EQ(m.getIndexRange(Entity::Face, Dir::X).getEnd(),
                  cells + MIndex::getUnitVector(0));
        EXPECT_EQ(m.getIndexRange(Entity::Face, Dir::Y).getBegin(), p0);
        EXPECT_EQ(m.getIndexRange(Entity::Face, Dir::Y).getEnd(),
                  cells + MIndex::getUnitVector(1));
    }
    { // physical domain
        EXPECT_EQ(m.getExtent()[0], domain[0]);
        EXPECT_EQ(m.getExtent()[1], domain[1]);
        EXPECT_EQ(m.getVolume(), (end - start).prod());
        EXPECT_EQ(m.getOrigin(), start);
        EXPECT_EQ(m.getGlobalOrigin(), origin);
        EXPECT_EQ(m.getRange().getBegin(), start);
        EXPECT_EQ(m.getRange().getEnd(), end);
    }
    { // mesh hull
        EXPECT_FALSE(m.isSubMesh());
    }
    { // MultiIndex
        const MIndex ref{2, 1};
        size_t iref = ref[0] + m.getSize(Entity::Cell)[0] * ref[1];
        EXPECT_EQ(m.getMultiIndex(iref, Entity::Cell), ref);
        iref = ref[0] + m.getSize(Entity::Node)[0] * ref[1];
        EXPECT_EQ(m.getMultiIndex(iref, Entity::Node), ref);
        for (size_t i = 0; i < Mesh::Dim; ++i) {
            iref = ref[0] + m.getSize(Entity::Face, i)[0] * ref[1];
            EXPECT_EQ(m.getMultiIndex(iref, Entity::Face, i), ref);
        }

        // global (global index origin is 0 for this case)
        iref = ref[0] + m.getSize(Entity::Face, Dir::Y)[0] * ref[1];
        EXPECT_EQ(m.getGlobalIndex(iref, Entity::Face, Dir::Y), ref);
    }

    { // coordinates
        const MIndex ref{2, 1};
        // cells
        size_t iref = ref[0] + m.getSize(Entity::Cell)[0] * ref[1];
        PointType cref = start + h * (PointType(ref) + 0.5);
        EXPECT_EQ(m.getCoords(iref, Entity::Cell), cref);
        cref += origin;
        EXPECT_EQ(m.getGlobalCoords(iref, Entity::Cell), cref);
        // nodes
        iref = ref[0] + m.getSize(Entity::Node)[0] * ref[1];
        cref = start + h * PointType(ref);
        EXPECT_EQ(m.getCoords(iref, Entity::Node), cref);
        cref += origin;
        EXPECT_EQ(m.getGlobalCoords(iref, Entity::Node), cref);
        // faces
        for (size_t i = 0; i < Mesh::Dim; ++i) {
            iref = ref[0] + m.getSize(Entity::Face, i)[0] * ref[1];
            cref = start + h * (PointType(ref) +
                                0.5 * (1.0 - PointType::getUnitVector(i)));
            EXPECT_EQ(m.getCoords(iref, Entity::Face, i), cref);
            cref += origin;
            EXPECT_EQ(m.getGlobalCoords(iref, Entity::Face, i), cref);
        }

        // indices
        for (const auto &c : m.getIterator(Entity::Cell)) {
            PointType cref = start + h * (PointType(c) + 0.5);
            EXPECT_EQ(m.getCoords(c, Entity::Cell), cref);
            cref += origin;
            EXPECT_EQ(m.getGlobalCoords(c, Entity::Cell), cref);
        }

        // iterator
        auto iter = m.getIterator(Entity::Cell);
        for (auto it = iter.begin(); it != iter.end(); ++it) {
            PointType cref = start + h * (PointType(*it) + 0.5);
            EXPECT_EQ(m.getCoords(it), cref);
            EXPECT_EQ(m.getCoords(it.getFlatIndex(), it.getEntity()), cref);
            EXPECT_EQ(m.getCoords(it.getMultiIndex(), it.getEntity()), cref);
            EXPECT_EQ(m.getCoords(*it, it.getEntity()), cref);
            cref += origin;
            EXPECT_EQ(m.getGlobalCoords(it), cref);
            EXPECT_EQ(m.getGlobalCoords(it.getFlatIndex(), it.getEntity()),
                      cref);
            EXPECT_EQ(m.getGlobalCoords(it.getMultiIndex(), it.getEntity()),
                      cref);
            EXPECT_EQ(m.getGlobalCoords(*it, it.getEntity()), cref);
        }
    }
    { // cell volume and size
        const RealType V = h.prod();
        auto iter = m.getIterator(Entity::Cell);
        for (auto it = iter.begin(); it != iter.end(); ++it) {
            EXPECT_EQ(m.getCellVolume(it.getFlatIndex()), V);
            EXPECT_EQ(m.getCellVolume(*it), V);
            EXPECT_EQ(m.getCellVolume(it), V);
            EXPECT_EQ(m.getCellSize(it.getFlatIndex()), h);
            EXPECT_EQ(m.getCellSize(*it), h);
            EXPECT_EQ(m.getCellSize(it), h);
        }
    }
    { // cell surfaces
        PointType surf_sum(0);
        const RealType cell_vol = h.prod();
        const PointType nx_p = PointType::getUnitVector(Dir::X);
        const PointType ny_p = PointType::getUnitVector(Dir::Y);
        const PointType nx_m = -nx_p;
        const PointType ny_m = -ny_p;
        auto iter = m.getIterator(Entity::Cell);
        for (auto it = iter.begin(); it != iter.end(); ++it) {
            // left and bottom face index
            const auto fi_00 = *it;
            // right face index
            const auto fi_10 = *it + MIndex::getUnitVector(Dir::X);
            // top face index
            const auto fi_11 = *it + MIndex::getUnitVector(Dir::Y);

            // surface vectors for cell it
            const PointType S00 = m.getSurface(fi_00, *it, Dir::X); // left
            const PointType S01 = m.getSurface(fi_00, *it, Dir::Y); // bottom
            const PointType S10 = m.getSurface(fi_10, *it, Dir::X); // right
            const PointType S11 = m.getSurface(fi_11, *it, Dir::Y); // top

            surf_sum += S00;
            surf_sum += S01;
            surf_sum += S10;
            surf_sum += S11;

            EXPECT_EQ(S00.unit(), nx_m);
            EXPECT_EQ(S01.unit(), ny_m);
            EXPECT_EQ(S10.unit(), nx_p);
            EXPECT_EQ(S11.unit(), ny_p);

            EXPECT_EQ(m.getSurfaceNormal(fi_00, *it, Dir::X), nx_m); // left
            EXPECT_EQ(m.getSurfaceNormal(fi_00, *it, Dir::Y), ny_m); // bottom
            EXPECT_EQ(m.getSurfaceNormal(fi_10, *it, Dir::X), nx_p); // right
            EXPECT_EQ(m.getSurfaceNormal(fi_11, *it, Dir::Y), ny_p); // top

            EXPECT_EQ(m.getSurfaceArea(fi_00, *it, Dir::X), cell_vol / h[0]);
            EXPECT_EQ(m.getSurfaceArea(fi_00, *it, Dir::Y), cell_vol / h[1]);
            EXPECT_EQ(m.getSurfaceArea(fi_10, *it, Dir::X), cell_vol / h[0]);
            EXPECT_EQ(m.getSurfaceArea(fi_11, *it, Dir::Y), cell_vol / h[1]);
        }
        EXPECT_EQ(surf_sum.norm(), 0);
    }
}

TEST(StructuredUniform, SubMesh)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<float, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;

    const PointType end(1);
    const MIndex cells(8);
    const PointType h = end / PointType(cells);
    Mesh m0(end, cells, MeshIntegrity::FullMesh);

    // boxes
    { // arbitrary (round up)
        const PointType sstart(0.2);
        const PointType send(0.6);
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1 + 1) * h);
    }
    { // matching cell boundaries (no rounding)
        const PointType sstart(0.25);
        const PointType send(0.5);
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1) * h);
    }
    // lower-dimensional extraction
    { // arbitrary interior sub-slice
        PointType sstart(0.2);
        PointType send(0.6);
        send[0] = sstart[0];
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1 + 1) * h);
    }
    { // arbitrary sub-slice on left boundary
        PointType sstart(0.2);
        PointType send(0.6);
        send[0] = sstart[0] = 0;
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1 + 1) * h);
    }
    { // arbitrary sub-slice in left boundary cell
        PointType sstart(0.2);
        PointType send(0.6);
        send[0] = sstart[0] = 0.5 * h[0];
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1 + 1) * h);
    }
    { // arbitrary sub-slice on right boundary
        PointType sstart(0.2);
        PointType send(0.6);
        send[0] = sstart[0] = 1;
        const auto m1 = m0.getSubMesh(sstart, send);
        MIndex i0(sstart / h);
        MIndex i1(send / h);
        i0[0] -= 1;
        i1 += 1;    // internal round up
        i1[0] -= 1; // not for boundary normal
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1) * h);
    }
    { // arbitrary point
        const PointType sstart(0.2);
        const PointType send(sstart);
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1 + 1) * h);
        EXPECT_EQ(m1->getIndexRange(Cubism::EntityType::Cell).size(), 1);
    }
    { // point on cell boundary
        const PointType sstart(0.25);
        const PointType send(sstart);
        const auto m1 = m0.getSubMesh(sstart, send);
        const MIndex i0(sstart / h);
        const MIndex i1(send / h);
        EXPECT_EQ(m1->getRange().getBegin(), PointType(i0) * h);
        EXPECT_EQ(m1->getRange().getEnd(), PointType(i1 + 1) * h);
        EXPECT_EQ(m1->getIndexRange(Cubism::EntityType::Cell).size(), 1);
    }
}
} // namespace
