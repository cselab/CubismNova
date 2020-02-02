// File       : HDFDriver.h
// Created    : Sat Jan 25 2020 10:05:58 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF read/write driver
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef HDFDRIVER_H_I1VZGUDV
#define HDFDRIVER_H_I1VZGUDV

#include "Cubism/Common.h"
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
    void write(const std::string &,
               const std::string &,
               const FileDataType *,
               const Mesh &,
               const Cubism::EntityType,
               const typename Mesh::IndexRangeType,
               const size_t,
               const typename Mesh::PointType,
               const double,
               const bool) const;

    void read(const std::string &,
              FileDataType *,
              const typename Mesh::IndexRangeType,
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
    typename Mesh::MultiIndex rank_index;
    void write(const std::string &,
               const std::string &,
               const FileDataType *,
               const Mesh &,
               const Cubism::EntityType,
               const typename Mesh::IndexRangeType,
               const size_t,
               const typename Mesh::PointType,
               const double,
               const bool) const;

    void read(const std::string &,
              FileDataType *,
              const typename Mesh::IndexRangeType,
              const size_t) const;
};

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* HDFDRIVER_H_I1VZGUDV */
