// File       : BlockMeshTest.cpp
// Created    : Sun Jan 05 2020 11:36:50 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Block mesh / sub-mesh tests
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/Alloc/AlignedBlockAllocator.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Core/Index.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"
#include <vector>

namespace
{
using namespace Cubism;

TEST(BlockMesh, Field)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;

    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshIntegrity = typename Mesh::MeshIntegrity;
    using PointType = typename Mesh::PointType;
    using Entity = typename Mesh::EntityType;
    using Range = typename Mesh::RangeType;

    const MIndex block_cells(8);
    const MIndex nblocks(2);

    // global mesh
    const PointType end(1);
    const MIndex cells = nblocks * block_cells;
    Mesh m(end, cells, MeshIntegrity::FullMesh);

    // custom field state
    struct MyFieldState {
        MIndex block_index;
        Mesh *mesh;
    };
    using CellField = Block::CellField<double, Mesh::Dim, MyFieldState>;
    using FC = Block::FieldContainer<CellField>;

    const PointType block_extent = m.getExtent() / PointType(nblocks);
    FC fields;
    std::vector<Mesh *> mfields;
    std::vector<IRange> vfaces(Mesh::Dim);
    for (int bz = 0; bz < nblocks[2]; ++bz) {
        for (int by = 0; by < nblocks[1]; ++by) {
            for (int bx = 0; bx < nblocks[0]; ++bx) {
                const MIndex bi{bx, by, bz};
                const PointType bstart =
                    m.getBegin() + PointType(bi) * block_extent;
                const PointType bend = bstart + block_extent;
                const MIndex cells = block_cells;
                MIndex nodes = cells;
                for (size_t i = 0; i < Mesh::Dim; ++i) {
                    MIndex faces(cells);
                    if (bi[i] == nblocks[i] - 1) {
                        ++nodes[i];
                        ++faces[i];
                    }
                    vfaces[i] = IRange(faces);
                }
                MyFieldState fs = {};
                mfields.push_back(new Mesh(m.getGlobalRange(),
                                           Range(bstart, bend),
                                           IRange(cells),
                                           IRange(nodes),
                                           vfaces,
                                           MeshIntegrity::SubMesh));
                fs.block_index = bi;
                fs.mesh = mfields.back();
                fields.pushBack(new CellField(IRange(cells), fs));
            }
        }
    }

    size_t sum_cells = 0;
    size_t sum_nodes = 0;
    size_t sum_faces_x = 0;
    size_t sum_faces_y = 0;
    size_t sum_faces_z = 0;
    for (const auto f : fields) {
        const MyFieldState &fs = f->getState();
        const MIndex bi = fs.block_index;
        const Mesh &fm = *fs.mesh;
        const PointType morigin = m.getBegin() + PointType(bi) * block_extent;
        const PointType gorigin = m.getGlobalBegin();
        EXPECT_EQ(fm.getBegin(), morigin);
        EXPECT_EQ(fm.getGlobalBegin(), gorigin);
        sum_cells += fm.size(Entity::Cell);
        sum_nodes += fm.size(Entity::Node);
        sum_faces_x += fm.size(Entity::Face, Dir::X);
        sum_faces_y += fm.size(Entity::Face, Dir::Y);
        sum_faces_z += fm.size(Entity::Face, Dir::Z);
    }
    EXPECT_EQ(m.size(Entity::Cell), sum_cells);
    EXPECT_EQ(m.size(Entity::Node), sum_nodes);
    EXPECT_EQ(m.size(Entity::Face, Dir::X), sum_faces_x);
    EXPECT_EQ(m.size(Entity::Face, Dir::Y), sum_faces_y);
    EXPECT_EQ(m.size(Entity::Face, Dir::Z), sum_faces_z);

    for (auto m : mfields) {
        delete m;
    }
}
} // namespace
