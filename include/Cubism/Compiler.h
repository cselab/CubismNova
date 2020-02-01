// File       : Compiler.h
// Created    : Mon Jan 27 2020 02:55:01 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Compiler specific pre-processor settings
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef COMPILER_H_ZNMPGTMC
#define COMPILER_H_ZNMPGTMC

// https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
#if defined(_MSC_VER)
#define DISABLE_WARNING_PUSH __pragma(warning(push))
#define DISABLE_WARNING_POP __pragma(warning(pop))
#define DISABLE_WARNING(warningNumber)                                         \
    __pragma(warning(disable : warningNumber))

#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER DISABLE_WARNING(4100)
#define DISABLE_WARNING_UNREFERENCED_FUNCTION DISABLE_WARNING(4505)

#elif defined(__GNUC__) || defined(__clang__)
#define DO_PRAGMA(X) _Pragma(#X)
#define DISABLE_WARNING_PUSH DO_PRAGMA(GCC diagnostic push)
#define DISABLE_WARNING_POP DO_PRAGMA(GCC diagnostic pop)
#define DISABLE_WARNING(warningName)                                           \
    DO_PRAGMA(GCC diagnostic ignored #warningName)

#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER                          \
    DISABLE_WARNING(-Wunused-parameter)
#define DISABLE_WARNING_UNREFERENCED_FUNCTION                                  \
    DISABLE_WARNING(-Wunused-function)

#else
#define DISABLE_WARNING_PUSH
#define DISABLE_WARNING_POP
#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
#define DISABLE_WARNING_UNREFERENCED_FUNCTION
#endif

// skip C++ related OpenMPI code
#define OMPI_SKIP_MPICXX

#endif /* COMPILER_H_ZNMPGTMC */
