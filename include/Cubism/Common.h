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

#ifndef CUBISM_DIMENSION
#define CUBISM_DIMENSION 3
#endif /* CUBISM_DIMENSION */

#ifndef CUBISM_ALIGNMENT
#define CUBISM_ALIGNMENT 32
#endif /* CUBISM_ALIGNMENT */

// Boolean compile switches for commandline:
// CUBISM_32BIT_INDEX
// CUBISM_OPTIMIZED_FIELD_OP

static_assert(CUBISM_DIMENSION > 0, "CUBISM_DIMENSION must be > 0");

NAMESPACE_BEGIN(Cubism)
#if 1 == CUBISM_DIMENSION
/// @brief Convenience direction indexing (C: C=direction)
enum class Dir { X = 0 };
/// @brief Convenience rank-1 tensor components indexing (C: C=component)
using Vector = Dir;
/// @brief Convenience rank-2 tensor components indexing (CD: C=component;
/// D=direction)
enum class Tensor { XX = 0 };
#elif 2 == CUBISM_DIMENSION
/// @brief Convenience direction indexing (C: C=direction)
enum class Dir { X = 0, Y };
/// @brief Convenience rank-1 tensor components indexing (C: C=component)
using Vector = Dir;
/// @brief Convenience rank-2 tensor components indexing (CD: C=component;
/// D=direction)
enum class Tensor { XX = 0, XY, YX, YY };
#elif 3 == CUBISM_DIMENSION
/// @brief Convenience direction indexing (C: C=direction)
enum class Dir { X = 0, Y, Z };
/// @brief Convenience rank-1 tensor components indexing (C: C=component)
using Vector = Dir;
/// @brief Convenience rank-2 tensor components indexing (CD: C=component;
/// D=direction)
enum class Tensor { XX = 0, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };
#endif /* 1 == CUBISM_DIMENSION */

/// @brief Cubism entity type descriptor
///        Cell: Cell entity, coordinates map to cell center (default)
///        Node: Node entity (vertices)
///        Face: Face entity that is spanned by nodes and is boundary of a cell.
///              Coordinates map the face center.
///        Undefined: No association
enum class EntityType { Cell = 0, Node, Face, Undefined };
NAMESPACE_END(Cubism)

// FIXME: [fabianw@mavt.ethz.ch; 2019-04-01] Deprecated, will be removed
#ifdef _USE_NUMA_
#warning _USE_NUMA_ is deprecated, use CUBISM_USE_NUMA instead.
#define CUBISM_USE_NUMA
#endif

#endif /* COMMON_H_C0GYQV59 */
