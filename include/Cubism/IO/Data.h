// File       : Data.h
// Created    : Wed Jun 09 2021 05:03:03 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Data low-level IO routines
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef DATA_H_7SBKSOVX
#define DATA_H_7SBKSOVX

#include "Cubism/Block/Data.h"
#include "Cubism/Common.h"
#include "Cubism/IO/HDFDriver.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER

// TODO: [fabianw@mavt.ethz.ch; 2021-03-24] Documentation in docs (rst files
// missing)

/**
 * @ingroup IO
 * @brief Convenience HDF5 data writer using a uniform mesh in [0, 1]
 * @tparam FileDataType HDF file data type
 * @tparam Data Data type
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``field``
 * @param data Input block data
 * @param time Current time
 * @param create_xdmf Flag for XDMF wrapper
 *
 * @rst
 * Write the data carried by ``data`` to an HDF5 container file.  This function
 * writes plain (low-level) data and assumes a uniform mesh in [0, 1].
 * @endrst
 */
template <typename FileDataType, typename Data>
void DataWriteUniformHDF(const std::string &fname,
                         const std::string &aname,
                         const Data &data,
                         const double time = 0.0,
                         const bool create_xdmf = true)
{
#ifdef CUBISM_USE_HDF
    using Mesh = Mesh::StructuredUniform<double, Data::IndexRangeType::Dim>;
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    Mesh mesh(MIndex(0),
              MIndex(1),
              data.getIndexRange().getExtent(),
              Mesh::MeshIntegrity::FullMesh);
    const IRange file_span = data.getIndexRange();
    const MIndex file_extent = file_span.getExtent();
    constexpr size_t NComp = 1; // SoA data
    if (create_xdmf) {
        std::printf(
            "DataWriteUniformHDF: Allocating %.1f kB file buffer (%s)\n",
            file_extent.prod() * NComp * sizeof(FileDataType) / 1024.,
            fname.c_str());
    }
    FileDataType *buf = new FileDataType[file_extent.prod() * NComp];
    for (auto i : file_span) { // copy data into buffer (FileDataType may not be
                               // equal to DataType)
        buf[file_span.getFlatIndex(i)] = static_cast<FileDataType>(data[i]);
    }
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.file_span = file_span;
    hdf_driver.write(
        fname,
        aname,
        buf,
        mesh,
        Cubism::EntityType::Cell, // all data treated as cell-centered
        NComp,
        time,
        create_xdmf);
    delete[] buf;
#else
    std::fprintf(
        stderr, "DataWriteUniformHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

DISABLE_WARNING_POP

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* DATA_H_7SBKSOVX */
