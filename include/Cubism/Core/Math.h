// File       : Math.h
// Created    : Thu Apr 25 2019 04:58:51 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Math related declarations
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef MATH_H_ODQVPXYD
#define MATH_H_ODQVPXYD

#include "Core/Common.h"

#include <cmath>
#include <cstdlib>

NAMESPACE_BEGIN(Cubism)

/// @brief Generic square root function.  std::sqrt() is overloaded for integral
///        types T and returns a double.  This causes implicit casts if the
///        return type of a function is DataType which may be an integral type.
///        This version forces a link-time error for any unsafe type T.
template <typename T>
inline T mySqrt(T v);

/// @brief Specialization of mySqrt for float type
template <>
inline float mySqrt(float v)
{
    return std::sqrt(v);
}

/// @brief Specialization of mySqrt for double type
template <>
inline double mySqrt(double v)
{
    return std::sqrt(v);
}

/// @brief Specialization of mySqrt for long double type
template <>
inline long double mySqrt(long double v)
{
    return std::sqrt(v);
}

/// @brief Generic abs function that bypasses unsigned integral types.
template <typename T>
inline T myAbs(T v)
{
    return std::abs(v);
}

/// @brief Specialization of myAbs for unsigned char type
template <>
inline unsigned char myAbs(unsigned char v)
{
    return v;
}

/// @brief Specialization of myAbs for unsigned short type
template <>
inline unsigned short myAbs(unsigned short v)
{
    return v;
}

/// @brief Specialization of myAbs for unsigned long type
template <>
inline unsigned long myAbs(unsigned long v)
{
    return v;
}

/// @brief Specialization of myAbs for unsigned long long type
template <>
inline unsigned long long myAbs(unsigned long long v)
{
    return v;
}

/// @brief Specialization of myAbs for unsigned type
template <>
inline unsigned myAbs(unsigned v)
{
    return v;
}

NAMESPACE_END(Cubism)

#endif /* MATH_H_ODQVPXYD */
