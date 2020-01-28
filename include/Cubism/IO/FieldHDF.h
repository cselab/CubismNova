// File       : FieldHDF.h
// Created    : Sat Jan 25 2020 01:11:02 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines for field types
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef FIELDHDF_H_SZ10D4OW
#define FIELDHDF_H_SZ10D4OW

#include "Cubism/Common.h"
#include "Cubism/Compiler.h"
#include "Cubism/IO/FieldAOS.h"
#include "Cubism/IO/HDFDriver.h"
#include <cstdio>
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
 * The ``Field`` type may be ``Cubism::FieldClass::Scalar`` or
 * ``Cubism::FieldClass::Tensor``.
 *
 * .. todo:: example for sub-space
 * @endrst
 */
template <typename FileDataType, typename Field, typename Mesh>
void FieldWriteHDF(const std::string &fname,
                   const std::string &aname,
                   const Field &field,
                   const Mesh &mesh,
                   const double time = 0,
                   const size_t face_dir = 0,
                   const bool create_xdmf = true)
{
#ifdef CUBISM_USE_HDF
    static_assert(Field::Class == Cubism::FieldClass::Scalar ||
                      Field::Class == Cubism::FieldClass::Tensor,
                  "FieldWriteHDF: Unsupported Cubism::FieldClass");
    using IRange = typename Mesh::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    if (Mesh::Dim > 1 && Mesh::Dim <= 3) {
        const IRange irange =
            mesh.getIndexRange(Field::BlockDataType::EntityType, face_dir);
        const MIndex iextent = irange.getExtent();
        constexpr size_t NComp = Field::NComponents;
        if (create_xdmf) {
            std::printf("FieldWriteHDF: Allocating %.1f kB file buffer (%s)\n",
                        iextent.prod() * NComp * sizeof(FileDataType) / 1024.,
                        fname.c_str());
        }
        FileDataType *buf = new FileDataType[iextent.prod() * NComp];
        Field2AOS(field, irange, buf);
        HDFDriver<FileDataType, typename Mesh::BaseMesh, Mesh::Class>
            hdf_driver;
        hdf_driver.write(fname,
                         aname,
                         buf,
                         mesh,
                         Field::BlockDataType::EntityType,
                         Field::Class,
                         NComp,
                         time,
                         face_dir,
                         create_xdmf);
        delete[] buf;
    } else {
        std::fprintf(stderr,
                     "FieldWriteHDF: Not supported for Mesh::Dim = %zu\n",
                     Mesh::Dim);
    }
#else
    std::fprintf(
        stderr, "FieldWriteHDF: HDF not supported (%s)\n", fname.c_str());
#endif /* CUBISM_USE_HDF */
}

// TODO: [fabianw@mavt.ethz.ch; 2020-01-24] Read

DISABLE_WARNING_POP

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* FIELDHDF_H_SZ10D4OW */
