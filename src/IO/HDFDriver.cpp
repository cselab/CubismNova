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

    (H5open() < 0) ? H5Eprint1(stderr) : 0;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file_id = H5Fcreate(
        (fname + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
    if (create_xdmf && Mesh::Dim > 1 && Mesh::Dim <= 3 &&
        entity != Cubism::EntityType::Face) {
        XDMFDriver<FileDataType, Class> xdmf;
        xdmf.write(fname, aname, mesh, entity, NComp, time);
    }

    // dataspace for dataset
    hid_t fspace_id = H5Screate_simple(HDFDim, dimsZYXC, NULL);
    hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);

    // chunked dataset
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    (H5Pset_chunk(plist_id, HDFDim, countZYXC) < 0) ? H5Eprint1(stderr) : 0;
    hid_t dataset_id = H5Dcreate(file_id,
                                 "data",
                                 getH5T<FileDataType>(),
                                 fspace_id,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
    (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;

    // define the hyperslab
    fspace_id = H5Dget_space(dataset_id);
    (H5Sselect_hyperslab(
         fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL) < 0)
        ? H5Eprint1(stderr)
        : 0;

    // create property list for collective write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    (H5Dwrite(dataset_id,
              getH5T<FileDataType>(),
              mspace_id,
              fspace_id,
              plist_id,
              buf) < 0)
        ? H5Eprint1(stderr)
        : 0;
    (H5Sclose(mspace_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Dclose(dataset_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Fclose(file_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5close() < 0) ? H5Eprint1(stderr) : 0;
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

    (H5open() < 0) ? H5Eprint1(stderr) : 0;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file_id = H5Fopen((fname + ".h5").c_str(), H5F_ACC_RDONLY, plist_id);
    (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;

    // dataspace for dataset
    hid_t dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
    hid_t fspace_id = H5Dget_space(dataset_id);
    hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);
    (H5Sselect_hyperslab(
         fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL) < 0)
        ? H5Eprint1(stderr)
        : 0;

    // read
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    (H5Dread(dataset_id,
             getH5T<FileDataType>(),
             mspace_id,
             fspace_id,
             plist_id,
             buf) < 0)
        ? H5Eprint1(stderr)
        : 0;
    (H5Sclose(mspace_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Dclose(dataset_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5Fclose(file_id) < 0) ? H5Eprint1(stderr) : 0;
    (H5close() < 0) ? H5Eprint1(stderr) : 0;
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
    const auto null = data_span.getNullSpace();
    const bool has_data = (null.size() < Mesh::Dim);
    int myrank, master_write;
    MPI_Comm comm_io;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_split(comm, has_data, myrank, &comm_io);
    MPI_Comm_rank(comm_io, &master_write);
    if (has_data) {
        const MIndex dbegin = data_span.getBegin();
        const MIndex fbegin = file_span.getBegin();
        const MIndex rextent = data_span.getExtent();
        const MIndex fextent = file_span.getExtent();
        constexpr size_t HDFDim = Mesh::Dim + 1;
        hsize_t offsetZYXC[HDFDim] = {};
        hsize_t countZYXC[HDFDim];
        hsize_t dimsZYXC[HDFDim];
        for (size_t i = 0; i < Mesh::Dim; ++i) {
            assert(dbegin[Mesh::Dim - 1 - i] >= fbegin[Mesh::Dim - 1 - i]);
            offsetZYXC[i] =
                dbegin[Mesh::Dim - 1 - i] - fbegin[Mesh::Dim - 1 - i];
            countZYXC[i] = rextent[Mesh::Dim - 1 - i];
            dimsZYXC[i] = fextent[Mesh::Dim - 1 - i];
        }
        offsetZYXC[HDFDim - 1] = 0;
        countZYXC[HDFDim - 1] = NComp;
        dimsZYXC[HDFDim - 1] = NComp;

        (H5open() < 0) ? H5Eprint1(stderr) : 0;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        (H5Pset_fapl_mpio(plist_id, comm_io, MPI_INFO_NULL) < 0)
            ? H5Eprint1(stderr)
            : 0;
        hid_t file_id = H5Fcreate(
            (fname + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
        if (master_write && create_xdmf && Mesh::Dim > 1 && Mesh::Dim <= 3 &&
            entity != Cubism::EntityType::Face) {
            XDMFDriver<FileDataType, Class> xdmf;
            xdmf.write(fname, aname, global_mesh, entity, NComp, time);
        }

        // dataspace for dataset
        hid_t fspace_id = H5Screate_simple(HDFDim, dimsZYXC, NULL);
        hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);

        // chunked dataset
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        (H5Pset_chunk(plist_id, HDFDim, countZYXC) < 0) ? H5Eprint1(stderr) : 0;
        hid_t dataset_id = H5Dcreate(file_id,
                                     "data",
                                     getH5T<FileDataType>(),
                                     fspace_id,
                                     H5P_DEFAULT,
                                     plist_id,
                                     H5P_DEFAULT);
        (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;

        // define the hyperslab
        fspace_id = H5Dget_space(dataset_id);
        (H5Sselect_hyperslab(
             fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL) < 0)
            ? H5Eprint1(stderr)
            : 0;

        // create property list for collective write
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        (H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) < 0)
            ? H5Eprint1(stderr)
            : 0;
        (H5Dwrite(dataset_id,
                  getH5T<FileDataType>(),
                  mspace_id,
                  fspace_id,
                  plist_id,
                  buf) < 0)
            ? H5Eprint1(stderr)
            : 0;
        (H5Dclose(dataset_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Sclose(mspace_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Fclose(file_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5close() < 0) ? H5Eprint1(stderr) : 0;
    }
    MPI_Comm_free(&comm_io);
}

template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
void HDFDriverMPI<FileDataType, Mesh, Class>::read(
    const std::string &fname,
    FileDataType *buf,
    const size_t NComp) const
{
    using MIndex = typename Mesh::MultiIndex;
    const auto null = data_span.getNullSpace();
    const bool has_data = (null.size() < Mesh::Dim);
    int myrank, master_write;
    MPI_Comm comm_io;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_split(comm, has_data, myrank, &comm_io);
    MPI_Comm_rank(comm_io, &master_write);
    if (has_data) {
        const MIndex dbegin = data_span.getBegin();
        const MIndex fbegin = file_span.getBegin();
        const MIndex rextent = data_span.getExtent();
        constexpr size_t HDFDim = Mesh::Dim + 1;
        hsize_t offsetZYXC[HDFDim] = {};
        hsize_t countZYXC[HDFDim];
        for (size_t i = 0; i < Mesh::Dim; ++i) {
            assert(dbegin[Mesh::Dim - 1 - i] >= fbegin[Mesh::Dim - 1 - i]);
            offsetZYXC[i] =
                dbegin[Mesh::Dim - 1 - i] - fbegin[Mesh::Dim - 1 - i];
            countZYXC[i] = rextent[Mesh::Dim - 1 - i];
        }
        offsetZYXC[HDFDim - 1] = 0;
        countZYXC[HDFDim - 1] = NComp;

        (H5open() < 0) ? H5Eprint1(stderr) : 0;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        (H5Pset_fapl_mpio(plist_id, comm_io, MPI_INFO_NULL) < 0)
            ? H5Eprint1(stderr)
            : 0;
        hid_t file_id =
            H5Fopen((fname + ".h5").c_str(), H5F_ACC_RDONLY, plist_id);
        (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;

        // dataspace for dataset
        hid_t dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
        hid_t fspace_id = H5Dget_space(dataset_id);
        hid_t mspace_id = H5Screate_simple(HDFDim, countZYXC, NULL);
        (H5Sselect_hyperslab(
             fspace_id, H5S_SELECT_SET, offsetZYXC, NULL, countZYXC, NULL) < 0)
            ? H5Eprint1(stderr)
            : 0;

        // collective read
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        (H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) < 0)
            ? H5Eprint1(stderr)
            : 0;
        (H5Dread(dataset_id,
                 getH5T<FileDataType>(),
                 mspace_id,
                 fspace_id,
                 plist_id,
                 buf) < 0)
            ? H5Eprint1(stderr)
            : 0;
        (H5Dclose(dataset_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Sclose(fspace_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Sclose(mspace_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Pclose(plist_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5Fclose(file_id) < 0) ? H5Eprint1(stderr) : 0;
        (H5close() < 0) ? H5Eprint1(stderr) : 0;
    }
    MPI_Comm_free(&comm_io);
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
                              Cubism::MeshClass::Uniform>;                     \
    template struct HDFDriverMPI<HDFType,                                      \
                                 Mesh::StructuredBase<float, Dim>,             \
                                 Cubism::MeshClass::Uniform>;                  \
    template struct HDFDriverMPI<HDFType,                                      \
                                 Mesh::StructuredBase<double, Dim>,            \
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
