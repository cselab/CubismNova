// File       : Version.h
// Created    : Wed Jan 15 2020 08:37:47 PM (+0100)
// Author     : Fabian Wermelinger
// Description: CubismNova build version strings
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef VERSION_H_BPFXCOEG
#define VERSION_H_BPFXCOEG

#include "Common.h"

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

/// @brief Cubism release version string
///
/// Command: git describe --abbrev=0
/// Example: v1.0.0, v1.0.0-rc1
extern const char *CubismVersion;

/// @brief Cubism HEAD at build time
///
/// Command: git describe --long --dirty --broken
/// Example: v1.1.8-5-g824e676-dirty
/// The string starts with the version number, followed by the number of commits
/// ahead of tagged version (5 in this case), followed by the SHA of HEAD,
/// followed by the state of the working tree (if any).
extern const char *CubismVersionHEAD;

/// @brief Cubism build branch
///
/// Command: git rev-parse --abbrev-ref HEAD
/// Example: master
extern const char *CubismBranch;

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* VERSION_H_BPFXCOEG */
