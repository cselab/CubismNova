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
template <typename FileDataType, typename Grid, typename Mesh>
void CartesianWriteHDF(const std::string &fname,
                       const std::string &aname,
                       const Grid &grid,
                       const Mesh &mesh,
                       const double time = 0,
                       const size_t face_dir = 0,
                       const bool create_xdmf = true)
{
#ifdef CUBISM_USE_HDF
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    constexpr typename Cubism::EntityType entity = Grid::EntityType;
    constexpr size_t NComp = Grid::NComponents;

    const IRange irange = mesh.getIndexRange(entity, face_dir);
    const MIndex iextent = irange.getExtent();
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
        Field2AOS(bf, irange, buf, face_dir);
    }
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.write(fname,
                     aname,
                     buf,
                     mesh,
                     entity,
                     NComp,
                     time,
                     face_dir,
                     create_xdmf);
    delete[] buf;
#else
    std::fprintf(
        stderr, "CartesianWriteHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

DISABLE_WARNING_POP

/**
 * @ingroup IO
 * @brief Write Cartesian grid data to HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Grid Grid type
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``grid``
 * @param grid Input grid
 * @param time Current time
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 * @param create_xdmf Flag for XDMF wrapper
 *
 * Convenience wrapper to dump a full grid to an HDF container file.
 */
template <typename FileDataType, typename Grid>
void CartesianWriteHDF(const std::string &fname,
                       const std::string &aname,
                       const Grid &grid,
                       const double time = 0,
                       const size_t face_dir = 0,
                       const bool create_xdmf = true)
{
    Cubism::IO::CartesianWriteHDF<FileDataType>(
        fname, aname, grid, grid.getMesh(), time, face_dir, create_xdmf);
}

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* CARTESIANHDF_H_GQ8CNN6W */
