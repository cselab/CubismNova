// File       : FieldHDF.h
// Created    : Sat Jan 25 2020 01:11:02 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines for field types
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef FIELDHDF_H_SZ10D4OW
#define FIELDHDF_H_SZ10D4OW

#include "Cubism/Common.h"
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

/**
 * @ingroup IO
 * @brief Write field data to HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Field Field type
 * @tparam Mesh Mesh type
 * @tparam Dir Special type that defines a cast to ``size_t``
 * @param fname Output full filename without file extension
 * @param aname Name of quantity in ``field``
 * @param field Input field
 * @param mesh Input mesh
 * @param time Current time
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 * @param create_xdmf Flag for XDMF wrapper
 *
 * @rst
 * Write the data carried by ``field`` to an HDF5 container file.  The data that
 * is written to the file is specified by the index space described in ``mesh``.
 *
 * .. note::
 *
 *    The variables ``field`` and ``mesh`` are only related by the index space
 *    they span.  If there is no common intersection, no data will be written to
 *    the file.  If the mesh spans a larger index space and full or partially
 *    contains the index space spanned by ``field`` then the index space spanned
 *    by ``mesh`` is clipped to that spanned by ``field``.
 *
 * Example:
 *
 * .. code-block:: cpp
 *
 *    #include "Cubism/Block/Field.h"
 *    #include "Cubism/Core/Index.h"
 *    #include "Cubism/Mesh/StructuredUniform.h"
 *    #include "Cubism/IO/FieldHDF.h"
 *
 *    using Mesh = Cubism::Mesh::StructuredUniform<double, 3>;
 *    using MeshIntegrity = typename Mesh::MeshIntegrity;
 *    using PointType = typename Mesh::PointType;
 *    using IRange = typename Mesh::IndexRangeType;
 *    using MIndex = typename Mesh::MultiIndex;
 *
 *    using CellField = Cubism::Block::CellField<double, Mesh::Dim>;
 *
 *    int main(void)
 *    {
 *        // mesh in [0,1]^3 with 16^3 cells
 *        const PointType end(1);
 *        const MIndex cells(16);
 *        Mesh m(end, cells, MeshIntegrity::FullMesh);
 *
 *        // cell field spans an index space [6,10)^3
 *        const IRange cell_range(6, 10);
 *        CellField cf(cell_range);
 *
 *        // write field to file (since the mesh fully contains cell_range, the
 *        // file contains 64 doubles at the coordinates defined in m using the
 *        // indices in cell_range
 *        Cubism::IO::FieldWriteHDF<double>(
 *            "filename", "attributename", cf, m, 0.0);
 *        return 0;
 *    }
 * @endrst
 */
template <typename FileDataType,
          typename Field,
          typename Mesh,
          typename Dir = size_t>
void FieldWriteHDF(const std::string &fname,
                   const std::string &aname,
                   const Field &field,
                   const Mesh &mesh,
                   const double time,
                   const Dir face_dir = 0,
                   const bool create_xdmf = true)
{
#ifdef CUBISM_USE_HDF
    static_assert(Field::Class == Cubism::FieldClass::Scalar ||
                      Field::Class == Cubism::FieldClass::Tensor ||
                      Field::Class == Cubism::FieldClass::FaceContainer,
                  "FieldReadHDF: Unsupported Cubism::FieldClass");
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    const size_t dface = static_cast<size_t>(face_dir);
    const auto clip =
        mesh.getSubMesh(field.getIndexRange(dface), Field::EntityType, dface);
    const IRange file_span = clip->getIndexRange(Field::EntityType, dface);
    const MIndex file_extent = file_span.getExtent();
    constexpr size_t NComp = Field::NComponents;
    if (create_xdmf) {
        std::printf("FieldWriteHDF: Allocating %.1f kB file buffer (%s)\n",
                    file_extent.prod() * NComp * sizeof(FileDataType) / 1024.,
                    fname.c_str());
    }
    FileDataType *buf = new FileDataType[file_extent.prod() * NComp];
    Field2AOS(field, file_span, buf, dface);
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.file_span = file_span;
    hdf_driver.write(fname,
                     aname,
                     buf,
                     *clip,
                     Field::EntityType,
                     NComp,
                     time,
                     create_xdmf);
    delete[] buf;
#else
    std::fprintf(
        stderr, "FieldWriteHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

/**
 * @ingroup IO
 * @brief Read field data from HDF file
 * @tparam FileDataType HDF file data type
 * @tparam Field Field type
 * @tparam Mesh Mesh type
 * @tparam Dir Special type that defines a cast to ``size_t``
 * @param fname Input full filename without file extension
 * @param field Field populated with file data
 * @param mesh Field (sub)mesh
 * @param face_dir Face direction (relevant for ``Cubism::EntityType::Face``)
 *
 * @rst
 * Read the data of an HDF5 container file into ``field``.  The data that is
 * read from the file is specified by the index space described in ``mesh``.
 * @endrst
 */
template <typename FileDataType,
          typename Field,
          typename Mesh,
          typename Dir = size_t>
void FieldReadHDF(const std::string &fname,
                  Field &field,
                  const Mesh &mesh,
                  const Dir face_dir = 0)
{
#ifdef CUBISM_USE_HDF
    static_assert(Field::Class == Cubism::FieldClass::Scalar ||
                      Field::Class == Cubism::FieldClass::Tensor ||
                      Field::Class == Cubism::FieldClass::FaceContainer,
                  "FieldReadHDF: Unsupported Cubism::FieldClass");
    {
        std::ifstream file(fname + ".h5");
        if (!file.good()) {
            throw std::runtime_error("FieldReadHDF: File '" + fname +
                                     "' does not exist");
        }
    }
    // XXX: [fabianw@mavt.ethz.ch; 2020-01-28] Instead of a mesh, an IndexRange
    // would be sufficient for the read operation; at the cost of an
    // inhomogeneous interface between FieldWriteHDF and FieldReadHDF.  I prefer
    // the mesh variant for this reason.
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    const size_t dface = static_cast<size_t>(face_dir);
    const auto clip =
        mesh.getSubMesh(field.getIndexRange(dface), Field::EntityType, dface);
    const IRange file_span = clip->getIndexRange(Field::EntityType, dface);
    const MIndex file_extent = file_span.getExtent();
    constexpr size_t NComp = Field::NComponents;
    FileDataType *buf = new FileDataType[file_extent.prod() * NComp];
    HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class> hdf_driver;
    hdf_driver.file_span = file_span;
    hdf_driver.read(fname, buf, NComp);
    AOS2Field(buf, file_span, field, dface);
    delete[] buf;
#else
    std::fprintf(
        stderr, "FieldReadHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

DISABLE_WARNING_POP

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* FIELDHDF_H_SZ10D4OW */
