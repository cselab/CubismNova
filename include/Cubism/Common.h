// File       : Common.h
// Created    : Mon Apr 01 2019 04:40:10 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Common declarations
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef COMMON_H_C0GYQV59
#define COMMON_H_C0GYQV59

#if !defined(NAMESPACE_BEGIN)
#define NAMESPACE_BEGIN(name)                                                  \
    namespace name                                                             \
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

static_assert(CUBISM_DIMENSION > 0, "CUBISM_DIMENSION must be > 0");

NAMESPACE_BEGIN(Cubism)
#if 1 == CUBISM_DIMENSION
/**
 * @brief Convenience direction indexing (C: C=direction)
 */
enum class Dir { X = 0 };
/**
 * @brief Convenience rank-1 tensor components indexing (C: C=component)
 */
using Vector = Dir;
/**
 * @brief Convenience rank-2 tensor components indexing (CD: C=component;
 * D=direction)
 */
enum class Tensor { XX = 0 };
#elif 2 == CUBISM_DIMENSION
/**
 * @brief Convenience direction indexing (C: C=direction)
 */
enum class Dir { X = 0, Y };
/**
 * @brief Convenience rank-1 tensor components indexing (C: C=component)
 */
using Vector = Dir;
/**
 * @brief Convenience rank-2 tensor components indexing (CD: C=component;
 * D=direction)
 */
enum class Tensor { XX = 0, XY, YX, YY };
#elif 3 == CUBISM_DIMENSION
/**
 * @brief Convenience direction indexing (C: C=direction)
 */
enum class Dir { X = 0, Y, Z };
/**
 * @brief Convenience rank-1 tensor components indexing (C: C=component)
 */
using Vector = Dir;
/**
 * @brief Convenience rank-2 tensor components indexing (CD: C=component;
 * D=direction)
 */
enum class Tensor { XX = 0, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };
#endif /* 1 == CUBISM_DIMENSION */

/**
 * @brief Cubism entity type descriptor
 *
 * @rst
 * Cell
 *    ``Cell`` entity, coordinates map to cell center (default)
 *
 * Node
 *    ``Node`` entity (vertices)
 *
 * Face
 *    ``Face`` entity that is spanned by nodes and is boundary of a cell
 *
 *    Coordinates map the face center
 *
 * Undefined
 *    No association
 * @endrst
 */
enum class EntityType { Cell = 0, Node, Face, Undefined };

/** @brief Field class identifier */
enum class FieldClass { Scalar = 0, Tensor, FaceContainer };

/** @brief Mesh class descriptor */
enum class MeshClass { Uniform = 0, Stretched };

/** @brief Mesh integrity type
 *
 * @rst
 * A ``FullMesh`` is one that describes the full domain local to a process.
 * A ``SubMesh`` is one that describes a sub-region of a ``FullMesh``.
 * @endrst */
enum class MeshIntegrity { FullMesh = 0, SubMesh };
NAMESPACE_END(Cubism)

#endif /* COMMON_H_C0GYQV59 */
