// File       : CartesianHDFTest.cpp
// Created    : Wed Jan 29 2020 10:53:13 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian grid HDF IO
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/IO/CartesianHDF.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"
#include <algorithm>

namespace
{
using namespace Cubism;

template <Cubism::EntityType Entity>
struct Initializer {
    template <typename Grid>
    void init(Grid &g)
    {
        typename Grid::DataType k = 0;
        for (auto f : g) {
            std::fill(f->begin(), f->end(), k++);
        }
    }
};

template <>
struct Initializer<Cubism::EntityType::Face> {
    template <typename Grid>
    void init(Grid &g)
    {
        typename Grid::DataType k = 0;
        for (auto ff : g) {      // face container
            for (auto f : *ff) { // face direction
                std::fill(f->begin(), f->end(), k++);
            }
        }
    }
};

TEST(IO, CartesianWriteHDF)
{
    {
        // 2D mesh
        using Mesh = Mesh::StructuredUniform<double, 2>;
        using MIndex = typename Mesh::MultiIndex;
        using PointType = typename Mesh::PointType;

        // grid blocks and cells per block
        const MIndex nblocks(3);
        const MIndex block_cells(8);

        using DataType = int;
        using Grid0 =
            Grid::Cartesian<DataType, Mesh, Cubism::EntityType::Cell, 0>;
        using Grid1 =
            Grid::Cartesian<DataType, Mesh, Cubism::EntityType::Node, 0>;
        using Grid2 =
            Grid::Cartesian<DataType, Mesh, Cubism::EntityType::Face, 0>;

        Grid0 grid(nblocks, block_cells);
        Initializer<Grid0::EntityType> ginit;
        ginit.init(grid);

        // all of the grid
        const auto &m = grid.getMesh();

        // sub-mesh selections
        //
        // sub-block (from 0.3 to 0.7)
        const auto blk = m.getSubMesh(PointType(0.3), PointType(0.7));
        // interior slice at x=0.3 (from 0.3 to 0.7)
        const auto ice = m.getSubMesh(PointType(0.3), PointType{0.3, 0.7});
        // larger sub-mesh request will default to same mesh as m
        const auto xxl = m.getSubMesh(PointType(-10.0), PointType(10.0));
        // slice at y=0 (left boundary)
        const auto lfy =
            m.getSubMesh(PointType{-10.0, 0.0}, PointType{10.0, 0.0});
        // slice at y=1 (right boundary)
        const auto rty =
            m.getSubMesh(PointType{-10.0, 1.0}, PointType{10.0, 1.0});
        // line at x=0.5 along y-direction
        const auto liz =
            m.getSubMesh(PointType{0.5, -10.0}, PointType{0.5, 10.0});

        const double time = 0;

        { // cell centered grid
            IO::CartesianWriteHDF<DataType>("cg2Dall", "basket", grid, time);
            IO::CartesianWriteHDF<DataType>(
                "cg2Dblk", "apples", grid, *blk, time);
            IO::CartesianWriteHDF<DataType>(
                "cg2Dice", "bananas", grid, *ice, time);
            IO::CartesianWriteHDF<DataType>(
                "cg2Dxxl", "peaches", grid, *xxl, time);
            IO::CartesianWriteHDF<DataType>(
                "cg2Dlfy", "oranges", grid, *lfy, time);
            IO::CartesianWriteHDF<DataType>(
                "cg2Drty", "kiwis", grid, *rty, time);
            IO::CartesianWriteHDF<DataType>(
                "cg2Dliz", "beers", grid, *liz, time);
        }
        { // node centered grid
            Grid1 grid(nblocks, block_cells);
            Initializer<Grid1::EntityType> ginit;
            ginit.init(grid);
            IO::CartesianWriteHDF<DataType>("ng2Dall", "basket", grid, time);
            IO::CartesianWriteHDF<DataType>(
                "ng2Dblk", "apples", grid, *blk, time);
            IO::CartesianWriteHDF<DataType>(
                "ng2Dice", "bananas", grid, *ice, time);
            IO::CartesianWriteHDF<DataType>(
                "ng2Dxxl", "peaches", grid, *xxl, time);
            IO::CartesianWriteHDF<DataType>(
                "ng2Dlfy", "oranges", grid, *lfy, time);
            IO::CartesianWriteHDF<DataType>(
                "ng2Drty", "kiwis", grid, *rty, time);
            IO::CartesianWriteHDF<DataType>(
                "ng2Dliz", "beers", grid, *liz, time);
        }
        { // face centered grid
            Grid2 grid(nblocks, block_cells);
            Initializer<Grid2::EntityType> ginit;
            ginit.init(grid);
            for (size_t d = 0; d < Grid2::Dim; ++d) {
                IO::CartesianWriteHDF<DataType>(
                    "fg2Dall" + std::to_string(d), "basket", grid, time, d);
                IO::CartesianWriteHDF<DataType>("fg2Dblk" + std::to_string(d),
                                                "apples",
                                                grid,
                                                *blk,
                                                time,
                                                d);
                IO::CartesianWriteHDF<DataType>("fg2Dice" + std::to_string(d),
                                                "bananas",
                                                grid,
                                                *ice,
                                                time,
                                                d);
                IO::CartesianWriteHDF<DataType>("fg2Dxxl" + std::to_string(d),
                                                "peaches",
                                                grid,
                                                *xxl,
                                                time,
                                                d);
                IO::CartesianWriteHDF<DataType>("fg2Dlfy" + std::to_string(d),
                                                "oranges",
                                                grid,
                                                *lfy,
                                                time,
                                                d);
                IO::CartesianWriteHDF<DataType>("fg2Drty" + std::to_string(d),
                                                "kiwis",
                                                grid,
                                                *rty,
                                                time,
                                                d);
                IO::CartesianWriteHDF<DataType>("fg2Dliz" + std::to_string(d),
                                                "beers",
                                                grid,
                                                *liz,
                                                time,
                                                d);
            }
        }
    }
    {
        // 3D mesh
        using Mesh = Mesh::StructuredUniform<double, 3>;
        using MIndex = typename Mesh::MultiIndex;
        using PointType = typename Mesh::PointType;

        // grid blocks and cells per block
        const MIndex nblocks(3);
        const MIndex block_cells(8);

        // Vector (rank-1) Cartesian node field
        using Grid3 = Grid::Cartesian<int, Mesh, Cubism::EntityType::Node, 1>;
        using DataType = typename Grid3::DataType;
        Grid3 grid(nblocks, block_cells);

        DataType k = 0;
        for (auto bf : grid) { // fill scalar block fields with data
            for (auto c : *bf) {
                std::fill(c->begin(), c->end(), k++);
            }
        }

        // all of the grid
        const auto &m = grid.getMesh();

        // sub-mesh selections
        const auto blk = m.getSubMesh(PointType(0.3), PointType(0.7));
        const auto ice = m.getSubMesh(PointType(0.3), PointType{0.3, 0.7, 0.7});
        const auto xxl = m.getSubMesh(PointType(-10.0), PointType(10.0));
        const auto lfy = m.getSubMesh(PointType{-10.0, 0.0, -10.0},
                                      PointType{10.0, 0.0, 10.0});
        const auto rty = m.getSubMesh(PointType{-10.0, 1.0, -10.0},
                                      PointType{10.0, 1.0, 10.0});
        const auto liz =
            m.getSubMesh(PointType{0.5, 0.5, -10.0}, PointType{0.5, 0.5, 10.0});

        const double time = 0;
        IO::CartesianWriteHDF<DataType>("ng3DVec", "basket", grid, time);
        IO::CartesianWriteHDF<DataType>(
            "ng3DVec_blk", "apples", grid, *blk, time);
        IO::CartesianWriteHDF<DataType>(
            "ng3DVec_ice", "bananas", grid, *ice, time);
        IO::CartesianWriteHDF<DataType>(
            "ng3DVec_xxl", "peaches", grid, *xxl, time);
        IO::CartesianWriteHDF<DataType>(
            "ng3DVec_lfy", "oranges", grid, *lfy, time);
        IO::CartesianWriteHDF<DataType>(
            "ng3DVec_rty", "kiwis", grid, *rty, time);
        IO::CartesianWriteHDF<DataType>(
            "ng3DVec_liz", "beers", grid, *liz, time);
    }
}

TEST(IO, CartesianWriteReadBackHDF)
{
    // 3D mesh
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using MIndex = typename Mesh::MultiIndex;

    // grid blocks and cells per block
    const MIndex nblocks(3);
    const MIndex block_cells(8);

    // Vector (rank-1) Cartesian cell field
    using Grid = Grid::Cartesian<int, Mesh, Cubism::EntityType::Cell, 0>;
    using DataType = typename Grid::DataType;
    Grid src(nblocks, block_cells);
    Initializer<Grid::EntityType> ginit;
    ginit.init(src);

    // write the field
    IO::CartesianWriteHDF<DataType>("framen", "ramen", src, 0);

    // read back soup
    Grid dst(nblocks, block_cells);
    IO::CartesianReadHDF<DataType>("framen", dst);

    int k = 0;
    for (const auto cf : dst) {
        for (const auto c : *cf) {
            EXPECT_EQ(c, k);
        }
        ++k;
    }
}
} // namespace
