// File       : Version.h
// Created    : Wed Jan 15 2020 08:37:47 PM (+0100)
// Author     : Fabian Wermelinger
// Description: CubismNova build version strings
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef VERSION_H_BPFXCOEG
#define VERSION_H_BPFXCOEG

#include "Cubism/Common.h"

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

/**
 * @ingroup Util
 * @brief Cubism release version string
 *
 * @rst
 * Command:
 *
 * .. code-block:: bash
 *
 *    git describe --abbrev=0
 *
 * Examples: ``v1.0.0``, ``v1.0.0-rc1``
 * @endrst
 */
extern const char *CubismVersion;

/**
 * @ingroup Util
 * @brief Cubism HEAD at build time
 *
 * @rst
 * Command:
 *
 * .. code-block:: bash
 *
 *    git describe --long --dirty --broken
 *
 * Example: ``v1.1.8-5-g824e676-dirty``
 *
 * The string starts with the version number, followed by the number of commits
 * ahead of tagged commit (``5`` in this case), followed by ``g`` and the short
 * hash (SHA-1)  of HEAD, followed by the state of the working tree (if any).
 * @endrst
 */
extern const char *CubismVersionHEAD;

/**
 * @ingroup Util
 * @brief Cubism build branch
 *
 * @rst
 * Command:
 *
 * .. code-block:: bash
 *
 *   git rev-parse --abbrev-ref HEAD
 *
 * Example: ``master``
 * @endrst
 */
extern const char *CubismBranch;

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* VERSION_H_BPFXCOEG */
