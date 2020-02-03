// File       : HDFDriver.h
// Created    : Sat Jan 25 2020 10:05:58 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF read/write driver
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef HDFDRIVER_H_I1VZGUDV
#define HDFDRIVER_H_I1VZGUDV

#include "Cubism/Common.h"
#include <mpi.h>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

/**
 * @brief HDF read/write interface
 * @tparam FileDataType File data taype
 * @tparam Mesh Mesh type
 * @tparam Class Mesh class
 * */
template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
struct HDFDriver {
    typename Mesh::IndexRangeType file_range;

    void write(const std::string &,
               const std::string &,
               const FileDataType *,
               const Mesh &,
               const Cubism::EntityType,
               const size_t,
               const double,
               const bool) const;

    void read(const std::string &,
              FileDataType *,
              const size_t) const;
};

/**
 * @brief HDF MPI read/write interface
 * @tparam FileDataType File data taype
 * @tparam Mesh Mesh type
 * @tparam Class Mesh class
 * */
template <typename FileDataType, typename Mesh, Cubism::MeshClass Class>
struct HDFDriverMPI {
    MPI_Comm comm;
    typename Mesh::IndexRangeType file_range;
    typename Mesh::IndexRangeType rank_range;

    void write(const std::string &,
               const std::string &,
               const FileDataType *,
               const Mesh &,
               const Cubism::EntityType,
               const size_t,
               const double,
               const bool) const;

    void read(const std::string &,
              FileDataType *,
              const size_t) const;
};

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* HDFDRIVER_H_I1VZGUDV */
