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
NAMESPACE_END(Cubism)

// FIXME: [fabianw@mavt.ethz.ch; 2019-04-01] Deprecated, will be removed
#ifdef _USE_NUMA_
#warning _USE_NUMA_ is deprecated, use CUBISM_USE_NUMA instead.
#define CUBISM_USE_NUMA
#endif

#ifdef _ON_FERMI_
#warning _ON_FERMI_ is deprecated, use CUBISM_ON_FERMI instead.
#define CUBISM_ON_FERMI
#endif

#ifdef _ALIGNBYTES_
#warning _ALIGNBYTES_ is deprecated, use CUBISM_ALIGNMENT instead.
#define CUBISM_ALIGNMENT _ALIGNBYTES_
#elif !defined(CUBISM_ALIGNMENT)
#define CUBISM_ALIGNMENT // If you get duplicate definition, put all Cubism
                         // includes after the main header include.
#endif

#endif /* COMMON_H_C0GYQV59 */
