// File       : HDFDriver.cpp
// Created    : Sat Jan 25 2020 09:09:07 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF IO routines implementation
// Copyright 2020 ETH Zurich. All Rights Reserved.

#include "Cubism/IO/HDFDriver.h"
#include "Cubism/Common.h"
#include "Cubism/IO/XDMFDriver.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

#ifdef CUBISM_USE_HDF
#include <hdf5.h>

template <typename T>
hid_t getH5T();

template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
void HDFDriver<FileDataType, Mesh, Class>::write(
    const std::string &fname,
    const std::string &aname,
    const FileDataType *buf,
    const Mesh &mesh,
    const Cubism::EntityType entity,
    const size_t NComp,
    const double time,
    const bool create_xdmf) const
{
    using MIndex = typename Mesh::MultiIndex;
    const MIndex fextent = file_span.getExtent(); // file range
    constexpr size_t HDFDim = Mesh::Dim + 1;
    hsize_t offsetZYXC[HDFDim] = {};
    hsize_t countZYXC[HDFDim];
    hsize_t dimsZYXC[HDFDim];
    for (size_t i = 0; i < Mesh::Dim; ++i) {
        countZYXC[i] = fextent[Mesh::Dim - 1 - i];
        dimsZYXC[i] = fextent[Mesh::Dim - 1 - i];
    }
    countZYXC[HDFDim - 1] = NComp;
    dimsZYXC[HDFDim - 1] = NComp;

    H5open();
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file_id =
        H5Fcreate((fname + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    (H5Pclose(fapl_id) < 0) ? H5Eprint1(stderr) : 0;
    if (create_xdmf && Mesh::Dim > 1 && Mesh::Dim <= 3 &&
        entity != Cubism::EntityType::Face) {
        XDMFDriver<FileDataType, Class> xdmf;
        xdmf.write(fname, aname, mesh, entity, NComp, time);
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

template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
void HDFDriver<FileDataType, Mesh, Class>::read(
    const std::string &fname,
    FileDataType *buf,
    const size_t NComp) const
{
    using MIndex = typename Mesh::MultiIndex;
    const MIndex fextent = file_span.getExtent(); // file range
    constexpr size_t HDFDim = Mesh::Dim + 1;
    hsize_t offsetZYXC[HDFDim] = {};
    hsize_t countZYXC[HDFDim];
    for (size_t i = 0; i < Mesh::Dim; ++i) {
        countZYXC[i] = fextent[Mesh::Dim - 1 - i];
    }
    countZYXC[HDFDim - 1] = NComp;

    H5open();
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file_id = H5Fopen((fname + ".h5").c_str(), H5F_ACC_RDONLY, fapl_id);
    (H5Pclose(fapl_id) < 0) ? H5Eprint1(stderr) : 0;
    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    hid_t dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
    hid_t fspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(
        fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL);
    hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);
    (H5Dread(dataset_id,
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

template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
void HDFDriverMPI<FileDataType, Mesh, Class>::write(
    const std::string &fname,
    const std::string &aname,
    const FileDataType *buf,
    const Mesh &global_mesh,
    const Cubism::EntityType entity,
    const size_t NComp,
    const double time,
    const bool create_xdmf) const
{
    using MIndex = typename Mesh::MultiIndex;
    const MIndex rbegin = rank_range.getBegin();
    const MIndex fbegin = file_range.getBegin();
    const MIndex rextent = rank_range.getExtent();
    const MIndex fextent = file_range.getExtent();
    constexpr size_t HDFDim = Mesh::Dim + 1;
    hsize_t offsetZYXC[HDFDim] = {};
    hsize_t countZYXC[HDFDim];
    hsize_t dimsZYXC[HDFDim];
    for (size_t i = 0; i < Mesh::Dim; ++i) {
        offsetZYXC[i] = rbegin[Mesh::Dim - 1 - i] - fbegin[Mesh::Dim - 1 - i];
        countZYXC[i] = rextent[Mesh::Dim - 1 - i];
        dimsZYXC[i] = fextent[Mesh::Dim - 1 - i];
    }
    offsetZYXC[HDFDim - 1] = 0;
    countZYXC[HDFDim - 1] = NComp;
    dimsZYXC[HDFDim - 1] = NComp;

    H5open();
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    (H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL) < 0) ? H5Eprint1(stderr)
                                                         : 0;
    hid_t file_id = H5Fopen((fname + ".h5").c_str(), H5F_ACC_RDWR, fapl_id);
    (H5Pclose(fapl_id) < 0) ? H5Eprint1(stderr) : 0;
    if (create_xdmf && Mesh::Dim > 1 && Mesh::Dim <= 3 &&
        entity != Cubism::EntityType::Face) {
        XDMFDriver<FileDataType, Class> xdmf;
        xdmf.write(fname, aname, global_mesh, entity, NComp, time);
    }
    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
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
    const auto null = rank_range.getNullSpace();
    if (null.size() == Mesh::Dim) {
        H5Sselect_none(fspace_id);
        H5Sselect_none(mspace_id);
    }
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

template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
void HDFDriverMPI<FileDataType, Mesh, Class>::read(
    const std::string &fname,
    FileDataType *buf,
    const size_t NComp) const
{
    using MIndex = typename Mesh::MultiIndex;
    const MIndex rbegin = rank_range.getBegin();
    const MIndex fbegin = file_range.getBegin();
    const MIndex rextent = rank_range.getExtent();
    const MIndex fextent = file_range.getExtent();
    constexpr size_t HDFDim = Mesh::Dim + 1;
    hsize_t offsetZYXC[HDFDim] = {};
    hsize_t countZYXC[HDFDim];
    hsize_t dimsZYXC[HDFDim];
    for (size_t i = 0; i < Mesh::Dim; ++i) {
        offsetZYXC[i] = rbegin[Mesh::Dim - 1 - i] - fbegin[Mesh::Dim - 1 - i];
        countZYXC[i] = rextent[Mesh::Dim - 1 - i];
        dimsZYXC[i] = fextent[Mesh::Dim - 1 - i];
    }
    offsetZYXC[HDFDim - 1] = 0;
    countZYXC[HDFDim - 1] = NComp;
    dimsZYXC[HDFDim - 1] = NComp;

    H5open();
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    (H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL) < 0) ? H5Eprint1(stderr)
                                                         : 0;
    hid_t file_id = H5Fopen((fname + ".h5").c_str(), H5F_ACC_RDONLY, fapl_id);
    (H5Pclose(fapl_id) < 0) ? H5Eprint1(stderr) : 0;
    hid_t dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
    hid_t fspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(
        fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL);
    hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);
    const auto null = rank_range.getNullSpace();
    if (null.size() == Mesh::Dim) {
        H5Sselect_none(fspace_id);
        H5Sselect_none(mspace_id);
    }
    (H5Dread(dataset_id,
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

// HDF5 type specializations
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

// explicit instantiations for CUBISM_DIMENSION
#define UNIFORM_HDFDRIVER(HDFType, Dim)                                        \
    template struct HDFDriver<HDFType,                                         \
                              Mesh::StructuredBase<float, Dim>,                \
                              Cubism::MeshClass::Uniform>;                     \
    template struct HDFDriver<HDFType,                                         \
                              Mesh::StructuredBase<double, Dim>,               \
                              Cubism::MeshClass::Uniform>;
UNIFORM_HDFDRIVER(float, 2)
UNIFORM_HDFDRIVER(double, 2)
UNIFORM_HDFDRIVER(size_t, 2)
UNIFORM_HDFDRIVER(int, 2)
UNIFORM_HDFDRIVER(char, 2)

UNIFORM_HDFDRIVER(float, 3)
UNIFORM_HDFDRIVER(double, 3)
UNIFORM_HDFDRIVER(size_t, 3)
UNIFORM_HDFDRIVER(int, 3)
UNIFORM_HDFDRIVER(char, 3)
#undef UNIFORM_HDFDRIVER

#endif /* CUBISM_USE_HDF */

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)
