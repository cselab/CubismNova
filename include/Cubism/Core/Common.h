// File       : Common.h
// Created    : Mon Apr 01 2019 04:40:10 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Common declarations
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef COMMON_H_C0GYQV59
#define COMMON_H_C0GYQV59

#if !defined(NAMESPACE_BEGIN)
#define NAMESPACE_BEGIN(name) \
    namespace name            \
    {
#endif

#if !defined(NAMESPACE_END)
#define NAMESPACE_END(name) }
#endif

NAMESPACE_BEGIN(Cubism)
/// @brief Constants for expressing direction
enum class Dir { X = 0, Y, Z, XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ, Any };

/// @brief Data layout constants used to describe a special block-data layout.
///        Undefined: no associated mapping (default)
///        Cell: data mapped to cell centers
///        Node: data mapped to cell nodes (vertices)
///        Face: data mapped to cell faces
enum class DataMapping { Undefined = 0, Cell, Node, Face };
NAMESPACE_END(Cubism)

#ifndef CUBISM_DIMENSION
#define CUBISM_DIMENSION 3
#endif /* CUBISM_DIMENSION */

#ifndef CUBISM_ALIGNMENT
#define CUBISM_ALIGNMENT 32
#endif /* CUBISM_ALIGNMENT */

// #ifndef CUBISM_32BIT_INDEX
// #define CUBISM_32BIT_INDEX
// #endif /* CUBISM_32BIT_INDEX */

// FIXME: [fabianw@mavt.ethz.ch; 2019-04-01] Deprecated, will be removed
#ifdef _USE_NUMA_
#warning _USE_NUMA_ is deprecated, use CUBISM_USE_NUMA instead.
#define CUBISM_USE_NUMA
#endif

#endif /* COMMON_H_C0GYQV59 */
