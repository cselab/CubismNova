// File       : FieldOperator.h
// Created    : Thu Jan 02 2020 03:43:07 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Block field math operators
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef FIELDOPERATOR_H_0XUYE2ZQ
#define FIELDOPERATOR_H_0XUYE2ZQ

#include <cassert>

#define OP_FIELD(A, B, C, N, OP)                                               \
    do {                                                                       \
        for (size_t i = 0; i < (N); ++i) {                                     \
            (*(C)++) = (*(A)++)OP(*(B)++);                                     \
        }                                                                      \
    } while (0)

#define OP_SCALAR(A, B, C, N, OP)                                              \
    do {                                                                       \
        for (size_t i = 0; i < (N); ++i) {                                     \
            (*(C)++) = (*(A)++)OP(B);                                          \
        }                                                                      \
    } while (0)

#define RCP(A, B, C, N)                                                        \
    do {                                                                       \
        for (size_t i = 0; i < (N); ++i) {                                     \
            (*(C)++) = (B) * (1 / (*(A)++));                                   \
        }                                                                      \
    } while (0)

#define CHECK_PTR(SRC)                                                         \
    do {                                                                       \
        if (sizeof(DataType) <= 8) {                                           \
            const DataType *check = (SRC);                                     \
            for (size_t i = 0; i < n; ++i) {                                   \
                assert((*check) > 0 || (*check) < 0);                          \
                ++check;                                                       \
            }                                                                  \
        }                                                                      \
    } while (0)

template <typename DataType>
void fieldAdd(const DataType *__restrict src0,
              const DataType *__restrict src1,
              DataType *__restrict dst,
              const size_t n)
{
    OP_FIELD(src0, src1, dst, n, +);
    return;
}

template <typename DataType>
void fieldSub(const DataType *__restrict src0,
              const DataType *__restrict src1,
              DataType *__restrict dst,
              const size_t n)
{
    OP_FIELD(src0, src1, dst, n, -);
    return;
}

template <typename DataType>
void fieldMul(const DataType *__restrict src0,
              const DataType *__restrict src1,
              DataType *__restrict dst,
              const size_t n)
{
    OP_FIELD(src0, src1, dst, n, *);
    return;
}

template <typename DataType>
void fieldDiv(const DataType *__restrict src0,
              const DataType *__restrict src1,
              DataType *__restrict dst,
              const size_t n)
{
#ifndef NDEBUG
    CHECK_PTR(src1);
#endif /* NDEBUG */
    OP_FIELD(src0, src1, dst, n, /);
    return;
}

template <typename DataType>
void fieldAdd(const DataType *__restrict src0,
              const DataType src1,
              DataType *__restrict dst,
              const size_t n)
{
    OP_SCALAR(src0, src1, dst, n, +);
    return;
}

template <typename DataType>
void fieldSub(const DataType *__restrict src0,
              const DataType src1,
              DataType *__restrict dst,
              const size_t n)
{
    OP_SCALAR(src0, src1, dst, n, -);
    return;
}

template <typename DataType>
void fieldMul(const DataType *__restrict src0,
              const DataType src1,
              DataType *__restrict dst,
              const size_t n)
{
    OP_SCALAR(src0, src1, dst, n, *);
    return;
}

template <typename DataType>
void fieldDiv(const DataType *__restrict src0,
              const DataType src1,
              DataType *__restrict dst,
              const size_t n)
{
#ifndef NDEBUG
    if (sizeof(DataType) <= 8) {
        assert(src1 > 0 || src1 < 0);
    }
#endif /* NDEBUG */
    OP_SCALAR(src0, src1, dst, n, /);
    return;
}

template <typename DataType>
void fieldRcp(const DataType *__restrict src0,
              const DataType src1,
              DataType *__restrict dst,
              const size_t n)
{
#ifndef NDEBUG
    CHECK_PTR(src0);
#endif /* NDEBUG */
    RCP(src0, src1, dst, n);
    return;
}

#undef OP_FIELD
#undef OP_SCALAR
#undef RCP
#undef CHECK_PTR

#ifdef CUBISM_OPTIMIZED_FIELD_OP
#define TARGET_NAME(OP, TYPE) field##OP##2_##TYPE
#define CALLER_NAME(OP) field##OP

#define OPTIMIZED_FIELD(OP, TYPE)                                              \
    extern "C" void TARGET_NAME(OP, TYPE)(                                     \
        const TYPE *, const TYPE *, TYPE *, const size_t);                     \
    template <>                                                                \
    void CALLER_NAME(OP)(const TYPE *__restrict src0,                          \
                         const TYPE *__restrict src1,                          \
                         TYPE *__restrict dst,                                 \
                         const size_t n)                                       \
    {                                                                          \
        TARGET_NAME(OP, TYPE)(src0, src1, dst, n);                             \
        return;                                                                \
    }

OPTIMIZED_FIELD(Add, double)
OPTIMIZED_FIELD(Sub, double)
OPTIMIZED_FIELD(Mul, double)
OPTIMIZED_FIELD(Div, double)
OPTIMIZED_FIELD(Add, float)
OPTIMIZED_FIELD(Sub, float)
OPTIMIZED_FIELD(Mul, float)
OPTIMIZED_FIELD(Div, float)
// XXX: [fabianw@mavt.ethz.ch; 2020-01-02] probably not
// OPTIMIZED_FIELD(Add, int)
// OPTIMIZED_FIELD(Sub, int)
// OPTIMIZED_FIELD(Mul, int)
// OPTIMIZED_FIELD(Div, int)

#undef OPTIMIZED_FIELD
#undef TARGET_NAME
#define TARGET_NAME(OP, TYPE) field##OP##1_##TYPE

#define OPTIMIZED_SCALAR(OP, TYPE)                                             \
    extern "C" void TARGET_NAME(OP, TYPE)(                                     \
        const TYPE *, const TYPE, TYPE *, const size_t);                       \
    template <>                                                                \
    void CALLER_NAME(OP)(const TYPE *__restrict src0,                          \
                         const TYPE src1,                                      \
                         TYPE *__restrict dst,                                 \
                         const size_t n)                                       \
    {                                                                          \
        TARGET_NAME(OP, TYPE)(src0, src1, dst, n);                             \
        return;                                                                \
    }

OPTIMIZED_SCALAR(Add, double)
OPTIMIZED_SCALAR(Sub, double)
OPTIMIZED_SCALAR(Mul, double)
OPTIMIZED_SCALAR(Div, double)
OPTIMIZED_SCALAR(Add, float)
OPTIMIZED_SCALAR(Sub, float)
OPTIMIZED_SCALAR(Mul, float)
OPTIMIZED_SCALAR(Div, float)
// XXX: [fabianw@mavt.ethz.ch; 2020-01-02] probably not
// OPTIMIZED_SCALAR(Add, int)
// OPTIMIZED_SCALAR(Sub, int)
// OPTIMIZED_SCALAR(Mul, int)
// OPTIMIZED_SCALAR(Div, int)

#undef OPTIMIZED_SCALAR
#undef CALLER_NAME
#undef TARGET_NAME
#define TARGET_NAME(TYPE) fieldRcp_##TYPE

#define OPTIMIZED_RCP(TYPE)                                                    \
    extern "C" void TARGET_NAME(TYPE)(                                         \
        const TYPE *, const TYPE, TYPE *, const size_t);                       \
    template <>                                                                \
    void fieldRcp(const TYPE *__restrict src0,                                 \
                  const TYPE src1,                                             \
                  TYPE *__restrict dst,                                        \
                  const size_t n)                                              \
    {                                                                          \
        TARGET_NAME(TYPE)(src0, src1, dst, n);                                 \
        return;                                                                \
    }

OPTIMIZED_RCP(double)
OPTIMIZED_RCP(float)
#undef OPTIMIZED_RCP
#undef TARGET_NAME
#endif /* CUBISM_OPTIMIZED_FIELD_OP */

#endif /* FIELDOPERATOR_H_0XUYE2ZQ */
