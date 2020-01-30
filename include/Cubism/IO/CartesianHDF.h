// File       : CartesianHDF.h
// Created    : Wed Jan 29 2020 10:25:26 AM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines for Cartesian grid types
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIANHDF_H_GQ8CNN6W
#define CARTESIANHDF_H_GQ8CNN6W

#include "Cubism/Common.h"
#include "Cubism/Compiler.h"
#include "Cubism/IO/FieldAOS.h"
#include "Cubism/IO/HDFDriver.h"
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER

// TODO: [fabianw@mavt.ethz.ch; 2020-01-29] FaceContainer interface
// /**
//  * @ingroup IO
//  * @brief Write Cartesian grid data to HDF file
//  * @tparam FileDataType HDF file data type
//  * @tparam Grid Grid type
//  * @tparam Mesh Mesh type
//  * @param fname Output full filename without file extension
//  * @param aname Name of quantity in ``grid``
//  * @param grid Input grid
//  * @param time Current time
//  * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
//  * @param create_xdmf Flag for XDMF wrapper
//  *
//  * @rst
//  * Write the data carried by ``grid`` to an HDF5 container file.  The data
//  that
//  * is written to the file is specified by the index space described in
//  ``mesh``.
//  * @endrst
//  */
// template <typename FileDataType, typename Grid, typename Mesh>
// void CartesianWriteHDF(const std::string &fname,
//                        const std::string &aname,
//                        const Grid &grid,
//                        const Mesh &mesh,
//                        const double time = 0,
//                        const size_t face_dir = 0,
//                        const bool create_xdmf = true)

/**
 * @ingroup IO
 * @brief Write Cartesian grid data to HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Grid Grid type
 * @tparam Mesh Mesh type
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``grid``
 * @param grid Input grid
 * @param time Current time
 * @param create_xdmf Flag for XDMF wrapper
 *
 * @rst
 * Write the data carried by ``grid`` to an HDF5 container file.  The data that
 * is written to the file is specified by the index space described in ``mesh``.
 * This function signature is used for scalar and tensor fields.
 * @endrst
 */
template <typename FileDataType, typename Grid, typename Mesh>
void CartesianWriteHDF(const std::string &fname,
                       const std::string &aname,
                       const Grid &grid,
                       const Mesh &mesh,
                       const double time = 0,
                       const bool create_xdmf = true)
{
#ifdef CUBISM_USE_HDF
    static_assert(Grid::BaseType::Class == Cubism::FieldClass::Scalar ||
                      Grid::BaseType::Class == Cubism::FieldClass::Tensor,
                  "CartesianWriteHDF: Unsupported Cubism::FieldClass");
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    constexpr typename Cubism::EntityType entity =
        Grid::BaseType::BlockDataType::EntityType;
    const IRange irange = mesh.getIndexRange(entity);
    const MIndex iextent = irange.getExtent();
    constexpr size_t NComp = Grid::BaseType::NComponents;
    if (create_xdmf) {
        std::printf("CartesianWriteHDF: Allocating %.1f MB file buffer (%s)\n",
                    iextent.prod() * NComp * sizeof(FileDataType) / 1024. /
                        1024.,
                    fname.c_str());
    }
    FileDataType *buf = new FileDataType[iextent.prod() * NComp];
#pragma omp parallel for
    for (size_t i = 0; i < grid.size(); ++i) {
        const auto &bf = grid[i]; // block field
        Field2AOS(bf, irange, buf);
    }
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.write(fname,
                     aname,
                     buf,
                     mesh,
                     entity,
                     Grid::BaseType::Class,
                     NComp,
                     time,
                     0,
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
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``grid``
 * @param grid Input grid
 * @param time Current time
 * @param create_xdmf Flag for XDMF wrapper
 *
 * Convenience wrapper to dump the full grid.
 */
template <typename FileDataType, typename Grid>
void CartesianWriteHDF(const std::string &fname,
                       const std::string &aname,
                       const Grid &grid,
                       const double time = 0,
                       const bool create_xdmf = true)
{
    Cubism::IO::CartesianWriteHDF<FileDataType>(
        fname, aname, grid, grid.getMesh(), time, create_xdmf);
}

DISABLE_WARNING_POP

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* CARTESIANHDF_H_GQ8CNN6W */
