// File       : Math.h
// Created    : Thu Apr 25 2019 04:58:51 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Math related declarations
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef MATH_H_ODQVPXYD
#define MATH_H_ODQVPXYD

#include "Common.h"

#include <cmath>
#include <cstdlib>

NAMESPACE_BEGIN(Cubism)

/**
 * @brief Generic square root function.
 * @tparam T Data type
 * @param v Argument
 *
 * @rst
 * Wrapper around ``std::sqrt()`` to force link-time errors for any unsafe types
 * ``T``.(e.g. integral types).  The template is specialized for:
 *
 * * ``float``
 * * ``double``
 * * ``long double``
 * @endrst
 */
template <typename T>
inline T mySqrt(T v);

template <>
inline float mySqrt(float v)
{
    return std::sqrt(v);
}

template <>
inline double mySqrt(double v)
{
    return std::sqrt(v);
}

template <>
inline long double mySqrt(long double v)
{
    return std::sqrt(v);
}

/**
 * @brief Generic abs function
 * @tparam T Data type
 * @param v Argument
 *
 * @rst
 * Wrapper around ``std::abs()`` to bypass unsigned integral types.  The
 * template is specialized for:
 *
 * * ``unsigned char``
 * * ``unsigned short``
 * * ``unsigned long``
 * * ``unsigned long long``
 * * ``unsigned``
 * @endrst
 */
template <typename T>
inline T myAbs(T v)
{
    return std::abs(v);
}

template <>
inline unsigned char myAbs(unsigned char v)
{
    return v;
}

template <>
inline unsigned short myAbs(unsigned short v)
{
    return v;
}

template <>
inline unsigned long myAbs(unsigned long v)
{
    return v;
}

template <>
inline unsigned long long myAbs(unsigned long long v)
{
    return v;
}

template <>
inline unsigned myAbs(unsigned v)
{
    return v;
}

NAMESPACE_END(Cubism)

#endif /* MATH_H_ODQVPXYD */
