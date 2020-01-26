// File       : HDFDriver.h
// Created    : Sat Jan 25 2020 10:05:58 PM (+0100)
// Author     : Fabian Wermelinger
// Description: HDF read/write driver
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef HDFDRIVER_H_I1VZGUDV
#define HDFDRIVER_H_I1VZGUDV

#include "Common.h"
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

/**
 * @brief HDF read/write interface
 * @tparam FileDataType File data taype
 * @tparam Mesh Mesh type
 * */
template <typename FileDataType, typename Mesh>
struct HDFDriver {
    void write(const std::string &,
               const std::string &,
               const FileDataType *,
               const Mesh &,
               const Cubism::EntityType,
               const Cubism::FieldClass,
               const size_t,
               const double,
               const size_t,
               const bool) const;
};

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* HDFDRIVER_H_I1VZGUDV */
