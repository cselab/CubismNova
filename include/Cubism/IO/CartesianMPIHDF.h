// File       : CartesianMPIHDF.h
// Created    : Wed Jan 29 2020 10:25:26 AM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines for Cartesian MPI grid types
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIANMPIHDF_H_S5X4YWDT
#define CARTESIANMPIHDF_H_S5X4YWDT

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
 * @brief Write Cartesian MPI grid data to HDF file
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
 * Write the data carried by the MPI ``grid`` to an HDF5 container file.  The
 * data that is written to the file is specified by the index space described in
 * ``mesh``.
 * @endrst
 */
template <typename FileDataType,
          typename Grid,
          typename Mesh,
          typename Dir = size_t>
void CartesianMPIWriteHDF(const std::string &fname,
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
                  "CartesianMPIWriteHDF: Unsupported Cubism::FieldClass");
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    constexpr typename Cubism::EntityType entity = Grid::EntityType;
    constexpr size_t NComp = Grid::NComponents;
    const size_t dface = static_cast<size_t>(face_dir);
    const auto &rmesh = grid.getMesh(); // rank local mesh
    const auto clip_global = // clip 'mesh' to the global grid mesh boundary
        mesh.getSubMesh(rmesh.getGlobalBegin(), rmesh.getGlobalEnd());
    const auto clip_rank = clip_global->getSubMesh(
        rmesh.getIndexRange(entity, dface), entity, dface);
    const IRange file_span = clip_global->getIndexRange(entity, dface);
    const IRange data_span = clip_rank->getIndexRange(entity, dface);
    const MIndex rank_extent = data_span.getExtent();
    if (create_xdmf && grid.isRoot()) {
        std::printf(
            "CartesianMPIWriteHDF: Allocating %.1f GB file buffer (%s)\n",
            file_span.getExtent().prod() * NComp * sizeof(FileDataType) /
                1024. / 1024. / 1024.,
            fname.c_str());
    }
    FileDataType *buf = new FileDataType[rank_extent.prod() * NComp];
#pragma omp parallel for
        for (size_t i = 0; i < grid.size(); ++i) {
        const auto &bf = grid[i]; // block field
        Field2AOS(bf, data_span, buf, dface);
    }
    HDFDriverMPI<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.comm = grid.getCartComm();
    hdf_driver.file_span = file_span;
    hdf_driver.data_span = data_span;
    hdf_driver.write(
        fname, aname, buf, *clip_global, entity, NComp, time, create_xdmf);
    delete[] buf;
#else
    std::fprintf(stderr,
                 "CartesianMPIWriteHDF: HDF not supported (%s)\n",
                 fname.c_str());
#endif /* CUBISM_USE_HDF */
}

/**
 * @ingroup IO
 * @brief Write Cartesian MPI grid data to HDF file
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
 * Convenience wrapper to dump a full MPI grid to an HDF container file.
 */
template <typename FileDataType, typename Grid, typename Dir = size_t>
void CartesianMPIWriteHDF(const std::string &fname,
                          const std::string &aname,
                          const Grid &grid,
                          const double time,
                          const Dir face_dir = 0,
                          const bool create_xdmf = true)
{
    Cubism::IO::CartesianMPIWriteHDF<FileDataType>(
        fname,
        aname,
        grid,
        grid.getGlobalMesh(),
        time,
        static_cast<size_t>(face_dir),
        create_xdmf);
}

/**
 * @ingroup IO
 * @brief Read Cartesian MPI grid data from HDF file
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
 * Read the data of an HDF5 container file into the MPI ``grid``.  The data that
 * is read from the file is specified by the index space described in ``mesh``.
 * @endrst
 */
template <typename FileDataType,
          typename Grid,
          typename Mesh,
          typename Dir = size_t>
void CartesianMPIReadHDF(const std::string &fname,
                         Grid &grid,
                         const Mesh &mesh,
                         const Dir face_dir = 0)
{
#ifdef CUBISM_USE_HDF
    static_assert(Grid::BaseType::Class == Cubism::FieldClass::Scalar ||
                      Grid::BaseType::Class == Cubism::FieldClass::Tensor ||
                      Grid::BaseType::Class ==
                          Cubism::FieldClass::FaceContainer,
                  "CartesianMPIReadHDF: Unsupported Cubism::FieldClass");
    {
        std::ifstream file(fname + ".h5");
        if (grid.isRoot() && !file.good()) {
            throw std::runtime_error("CartesianMPIReadHDF: File '" + fname +
                                     "' does not exist");
        }
    }
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    constexpr typename Cubism::EntityType entity = Grid::EntityType;
    constexpr size_t NComp = Grid::NComponents;
    const size_t dface = static_cast<size_t>(face_dir);
    const auto &rmesh = grid.getMesh(); // rank local mesh
    const auto clip_global = // clip 'mesh' to the global grid mesh boundary
        mesh.getSubMesh(rmesh.getGlobalBegin(), rmesh.getGlobalEnd());
    const auto clip_rank = clip_global->getSubMesh(
        rmesh.getIndexRange(entity, dface), entity, dface);
    const IRange file_span = clip_global->getIndexRange(entity, dface);
    const IRange data_span = clip_rank->getIndexRange(entity, dface);
    const MIndex rank_extent = data_span.getExtent();
    FileDataType *buf = new FileDataType[rank_extent.prod() * NComp];
    HDFDriverMPI<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.comm = grid.getCartComm();
    hdf_driver.file_span = file_span;
    hdf_driver.data_span = data_span;
    hdf_driver.read(fname, buf, NComp);
#pragma omp parallel for
    for (size_t i = 0; i < grid.size(); ++i) {
        auto &bf = grid[i]; // block field
        AOS2Field(buf, data_span, bf, dface);
    }
    delete[] buf;
#else
    std::fprintf(
        stderr, "CartesianMPIReadHDF: HDF not supported (%s)\n", fname.c_str());
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
 * Convenience wrapper to read a full MPI grid from an HDF container file.
 */
template <typename FileDataType, typename Grid, typename Dir = size_t>
void CartesianMPIReadHDF(const std::string &fname,
                         Grid &grid,
                         const Dir face_dir = 0)
{
    Cubism::IO::CartesianMPIReadHDF<FileDataType>(
        fname, grid, grid.getGlobalMesh(), static_cast<size_t>(face_dir));
}

DISABLE_WARNING_POP

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* CARTESIANMPIHDF_H_S5X4YWDT */
