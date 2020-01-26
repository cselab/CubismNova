// File       : HDFDriver.cpp
// Created    : Sat Jan 25 2020 09:09:07 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines implementation
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "Common.h"
#include "IO/XDMFDriver.h"
#include "Mesh/StructuredUniform.h"
#include <cstddef>
#include <cstdio>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

#ifdef CUBISM_USE_HDF
#include <hdf5.h>

template <typename T>
hid_t getH5T();
#endif /* CUBISM_USE_HDF */

template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
struct HDFDriver {
    void write(const std::string &fname,
               const std::string &aname,
               const FileDataType *buf,
               const Mesh &mesh,
               const Cubism::EntityType entity,
               const Cubism::FieldClass fclass,
               const size_t NComp,
               const double time,
               const size_t face_dir,
               const bool create_xdmf) const
    {
#ifdef CUBISM_USE_HDF
        using IRange = typename Mesh::IndexRangeType;
        using MIndex = typename IRange::MultiIndex;
        const IRange irange = mesh.getIndexRange(entity, face_dir);
        const MIndex iextent = irange.getExtent();
        constexpr size_t HDFDim = Mesh::Dim + 1;
        hsize_t offsetZYXC[HDFDim] = {};
        hsize_t countZYXC[HDFDim];
        hsize_t dimsZYXC[HDFDim];
        for (size_t i = 0; i < Mesh::Dim; ++i) {
            countZYXC[i] = iextent[Mesh::Dim - 1 - i];
            dimsZYXC[i] = iextent[Mesh::Dim - 1 - i];
        }
        countZYXC[HDFDim - 1] = NComp;
        dimsZYXC[HDFDim - 1] = NComp;

        H5open();
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        hid_t file_id = H5Fcreate(
            (fname + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        (H5Pclose(fapl_id) < 0) ? H5Eprint1(stderr) : 0;
        if (create_xdmf && entity != Cubism::EntityType::Face) {
            // TODO: [fabianw@mavt.ethz.ch; 2020-01-25] face storage location is
            // currently not supported by ParaView
            XDMFDriver<FileDataType, Class> xdmf;
            xdmf.write(fname, aname, mesh, fclass, entity, NComp, time);
        }
        fapl_id = H5Pcreate(H5P_DATASET_XFER);
        hid_t fspace_id = H5Screate_simple(HDFDim, dimsZYXC, NULL);
        hid_t dataset_id = H5Dcreate(file_id,
                                     "data",
                                     getH5T<FileDataType>(),
                                     fspace_id,
                                     H5P_DEFAULT,
                                     H5P_DEFAULT,
                                     H5P_DEFAULT);
        fspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(
            fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL);
        hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);
        (H5Dwrite(dataset_id,
                  getH5T<FileDataType>(),
                  mspace_id,
                  fspace_id,
                  fapl_id,
                  buf) < 0)
            ? H5Eprint1(stderr)
            : 0;
        (H5Sclose(mspace_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Dclose(dataset_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Pclose(fapl_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Fclose(file_id) < 0) ? H5Eprint1(stderr) : 0;
        H5close();
    }
#endif /* CUBISM_USE_HDF */
};

// HDF5 type specializations
#ifdef CUBISM_USE_HDF
template <>
inline hid_t getH5T<float>()
{
    return H5T_NATIVE_FLOAT;
}
template <>
inline hid_t getH5T<double>()
{
    return H5T_NATIVE_DOUBLE;
}
template <>
inline hid_t getH5T<size_t>()
{
    return H5T_NATIVE_UINT64;
}
template <>
inline hid_t getH5T<int>()
{
    return H5T_NATIVE_INT;
}
template <>
inline hid_t getH5T<char>()
{
    return H5T_NATIVE_CHAR;
}
#endif /* CUBISM_USE_HDF */

// explicit instantiations for CUBISM_DIMENSION
#define UNIFORM_HDFDRIVER(HDFType, MeshReal)                                   \
    template struct HDFDriver<                                                 \
        HDFType,                                                               \
        Mesh::StructuredBase<MeshReal, CUBISM_DIMENSION>,                      \
        Cubism::MeshClass::Uniform>
UNIFORM_HDFDRIVER(float, float);
UNIFORM_HDFDRIVER(float, double);
UNIFORM_HDFDRIVER(double, float);
UNIFORM_HDFDRIVER(double, double);
UNIFORM_HDFDRIVER(size_t, float);
UNIFORM_HDFDRIVER(size_t, double);
UNIFORM_HDFDRIVER(int, float);
UNIFORM_HDFDRIVER(int, double);
UNIFORM_HDFDRIVER(char, float);
UNIFORM_HDFDRIVER(char, double);
#undef UNIFORM_HDFDRIVER

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)