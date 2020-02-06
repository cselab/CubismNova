// File       : CartesianHDF.h
// Created    : Wed Jan 29 2020 10:25:26 AM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines for Cartesian grid types
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIANHDF_H_GQ8CNN6W
#define CARTESIANHDF_H_GQ8CNN6W

#include "Cubism/Common.h"
#include "Cubism/IO/FieldAOS.h"
#include "Cubism/IO/HDFDriver.h"
#include <cstdio>
#include <fstream>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER

/**
 * @ingroup IO
 * @brief Write Cartesian grid data to HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Grid Grid type
 * @tparam Mesh Mesh type
 * @tparam Dir Special type that defines a cast to ``size_t``
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``grid``
 * @param grid Input grid
 * @param mesh Input mesh corresponding to the extracted data
 * @param time Current time
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 * @param create_xdmf Flag for XDMF wrapper
 *
 * @rst
 * Write the data carried by ``grid`` to an HDF5 container file.  The data that
 * is written to the file is specified by the index space described in ``mesh``.
 * @endrst
 */
template <typename FileDataType,
          typename Grid,
          typename Mesh,
          typename Dir = size_t>
void CartesianWriteHDF(const std::string &fname,
                       const std::string &aname,
                       const Grid &grid,
                       const Mesh &mesh,
                       const double time,
                       const Dir face_dir = 0,
                       const bool create_xdmf = true)
{
#ifdef CUBISM_USE_HDF
    static_assert(Grid::BaseType::Class == Cubism::FieldClass::Scalar ||
                      Grid::BaseType::Class == Cubism::FieldClass::Tensor ||
                      Grid::BaseType::Class ==
                          Cubism::FieldClass::FaceContainer,
                  "CartesianWriteHDF: Unsupported Cubism::FieldClass");
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    constexpr typename Cubism::EntityType entity = Grid::EntityType;
    constexpr size_t NComp = Grid::NComponents;
    const size_t dface = static_cast<size_t>(face_dir);
    const auto clip = mesh.getSubMesh(
        grid.getMesh().getIndexRange(entity, dface), entity, dface);
    const IRange file_span = clip->getIndexRange(entity, dface);
    const MIndex file_extent = file_span.getExtent();
    if (create_xdmf) {
        std::printf("CartesianWriteHDF: Allocating %.1f MB file buffer (%s)\n",
                    file_extent.prod() * NComp * sizeof(FileDataType) / 1024. /
                        1024.,
                    fname.c_str());
    }
    FileDataType *buf = new FileDataType[file_extent.prod() * NComp];
#pragma omp parallel for
    for (size_t i = 0; i < grid.size(); ++i) {
        const auto &bf = grid[i]; // block field
        Field2AOS(bf, file_span, buf, dface);
    }
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.file_span = file_span;
    hdf_driver.write(fname,
                     aname,
                     buf,
                     *clip,
                     entity,
                     NComp,
                     time,
                     create_xdmf);
    delete[] buf;
#else
    std::fprintf(
        stderr, "CartesianWriteHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

/**
 * @ingroup IO
 * @brief Write Cartesian grid data to HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Grid Grid type
 * @tparam Dir Special type that defines a cast to ``size_t``
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``grid``
 * @param grid Input grid
 * @param time Current time
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 * @param create_xdmf Flag for XDMF wrapper
 *
 * Convenience wrapper to dump a full grid to an HDF container file.
 */
template <typename FileDataType, typename Grid, typename Dir = size_t>
void CartesianWriteHDF(const std::string &fname,
                       const std::string &aname,
                       const Grid &grid,
                       const double time,
                       const Dir face_dir = 0,
                       const bool create_xdmf = true)
{
    Cubism::IO::CartesianWriteHDF<FileDataType>(fname,
                                                aname,
                                                grid,
                                                grid.getGlobalMesh(),
                                                time,
                                                static_cast<size_t>(face_dir),
                                                create_xdmf);
}

/**
 * @ingroup IO
 * @brief Read Cartesian grid data from HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Grid Grid type
 * @tparam Mesh Mesh type
 * @tparam Dir Special type that defines a cast to ``size_t``
 * @param fname Input full filename without file extension
 * @param grid Grid populated with file data
 * @param mesh Grid (sub)mesh
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 *
 * @rst
 * Read the data of an HDF5 container file into ``grid``.  The data that is
 * read from the file is specified by the index space described in ``mesh``.
 * @endrst
 */
template <typename FileDataType,
          typename Grid,
          typename Mesh,
          typename Dir = size_t>
void CartesianReadHDF(const std::string &fname,
                      Grid &grid,
                      const Mesh &mesh,
                      const Dir face_dir = 0)
{
#ifdef CUBISM_USE_HDF
    static_assert(Grid::BaseType::Class == Cubism::FieldClass::Scalar ||
                      Grid::BaseType::Class == Cubism::FieldClass::Tensor ||
                      Grid::BaseType::Class ==
                          Cubism::FieldClass::FaceContainer,
                  "CartesianReadHDF: Unsupported Cubism::FieldClass");
    {
        std::ifstream file(fname + ".h5");
        if (!file.good()) {
            throw std::runtime_error("CartesianReadHDF: File '" + fname +
                                     "' does not exist");
        }
    }
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    constexpr typename Cubism::EntityType entity = Grid::EntityType;
    constexpr size_t NComp = Grid::NComponents;
    const size_t dface = static_cast<size_t>(face_dir);
    const auto clip = mesh.getSubMesh(
        grid.getMesh().getIndexRange(entity, dface), entity, dface);
    const IRange file_span = clip->getIndexRange(entity, dface);
    const MIndex file_extent = file_span.getExtent();
    FileDataType *buf = new FileDataType[file_extent.prod() * NComp];
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.file_span = file_span;
    hdf_driver.read(fname, buf, NComp);
#pragma omp parallel for
    for (size_t i = 0; i < grid.size(); ++i) {
        auto &bf = grid[i]; // block field
        AOS2Field(buf, file_span, bf, dface);
    }
    delete[] buf;
#else
    std::fprintf(
        stderr, "CartesianReadHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

/**
 * @ingroup IO
 * @brief Read Cartesian grid data from HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Grid Grid type
 * @tparam Dir Special type that defines a cast to ``size_t``
 * @param fname Input full filename without file extension
 * @param grid Grid populated with file data
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 *
 * Convenience wrapper to read a full grid from an HDF container file.
 */
template <typename FileDataType, typename Grid, typename Dir = size_t>
void CartesianReadHDF(const std::string &fname,
                      Grid &grid,
                      const Dir face_dir = 0)
{
    Cubism::IO::CartesianReadHDF<FileDataType>(
        fname, grid, grid.getGlobalMesh(), static_cast<size_t>(face_dir));
}

DISABLE_WARNING_POP

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* CARTESIANHDF_H_GQ8CNN6W */
