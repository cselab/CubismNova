// File       : BlockMeshTest.cpp
// Created    : Sun Jan 05 2020 11:36:50 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Block mesh / sub-mesh tests
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Block/Field.h"
#include "Core/Index.h"
#include "gtest/gtest.h"

#include "Alloc/AlignedBlockAllocator.h"
#include "Mesh/StructuredUniform.h"

#include <vector>

namespace
{
using namespace Cubism;

template <typename T, template <typename> class TAlloc, size_t DIM>
using CellData = Block::Data<T, EntityType::Cell, DIM, TAlloc<T>>;

TEST(BlockMesh, Field)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;

    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using MeshHull = typename Mesh::MeshHull;
    using PointType = typename Mesh::PointType;
    using Entity = typename Mesh::EntityType;
    using Range = typename Mesh::RangeType;

    const MIndex block_cells(8);
    const MIndex nblocks(2);

    // global mesh
    const PointType end(1);
    const MIndex cells = nblocks * block_cells;
    Mesh m(end, cells, MeshHull::FullMesh);

    // custom field state (rank and comp are required)
    struct MyFieldState {
        size_t rank;
        size_t comp;
        MIndex idx;
        Mesh *mesh;
    };
    using CellField =
        Block::Field<CellData<double, AlignedBlockAllocator, Mesh::Dim>,
                     MyFieldState>;
    using FC = Block::FieldContainer<CellField>;

    const PointType block_extent = m.getExtent() / PointType(nblocks);
    FC fields;
    std::vector<Mesh *> mfields;
    std::vector<IRange> vfaces(Mesh::Dim);
    for (int bz = 0; bz < nblocks[2]; ++bz) {
        for (int by = 0; by < nblocks[1]; ++by) {
            for (int bx = 0; bx < nblocks[0]; ++bx) {
                const PointType gorigin = m.getGlobalOrigin();
                const MIndex bi{bx, by, bz};
                const PointType bstart =
                    m.getOrigin() + PointType(bi) * block_extent;
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
                mfields.push_back(new Mesh(gorigin,
                                           Range(bstart, bend),
                                           IRange(cells),
                                           IRange(nodes),
                                           vfaces,
                                           MeshHull::SubMesh));
                fs.idx = bi;
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
        const MIndex bi = fs.idx;
        const Mesh &fm = *fs.mesh;
        const PointType morigin = m.getOrigin() + PointType(bi) * block_extent;
        const PointType gorigin = m.getGlobalOrigin();
        EXPECT_EQ(fm.getOrigin(), morigin);
        EXPECT_EQ(fm.getGlobalOrigin(), gorigin);
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
