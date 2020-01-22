// File       : FieldAOS.h
// Created    : Tue Jan 21 2020 02:10:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Copy data from a block field into an array of structures buffer
//              (AoS)
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef FIELDAOS_H_DOY3MTA8
#define FIELDAOS_H_DOY3MTA8

#include "Block/Field.h"
#include "Common.h"
#include <cstddef>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)
/**
 * @addtogroup IO
 * @{
 */

/**
 * @brief Copy scalar field data into AoS buffer
 * @tparam U Field block data type
 * @tparam V Field state type
 * @param f Source scalar field
 * @param r Index space for the copy (describes the memory region of buf)
 * @param buf Destination buffer
 *
 * Copy the data from a structure of arrays (SoA) field into an array of
 * structures (AoS) buffer for I/O operation.  This is a low-level function
 * which can be used in high-level I/O interfaces.
 */
template <typename U, typename V>
void FieldAOS(const Block::Field<U, V> &f,
              const typename Block::Field<U, V>::IndexRangeType &r,
              typename Block::Field<U, V>::DataType *buf)
{
    using IRange = typename Block::Field<U, V>::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    // get index space intersection.  The index range r must be relative to the
    // memory region allocated in buf.
    const IRange indices = r.getIntersection(f.getIndexRange());

    // local block offset
    const MIndex base = indices.getBegin();
    const MIndex offset = base - f.getIndexRange().getBegin();

    // iterate over local (sub) index space and copy into buf
    for (const auto &p : indices) {
        const size_t i = r.getFlatIndexFromGlobal(base + p);
        buf[i] = f[offset + p];
    }
}

/**
 * @brief Copy scalar field data into AoS buffer
 * @tparam U Field type
 * @tparam R Tensor rank
 * @param tf Source tensor field
 * @param r Index space for the copy (describes the memory region of buf)
 * @param buf Destination buffer
 *
 * Copy the data from a structure of arrays (SoA) field into an array of
 * structures (AoS) buffer for I/O operation.  This is a low-level function
 * which can be used in high-level I/O interfaces.
 */
template <typename U, size_t R>
void FieldAOS(const Block::TensorField<U, R> &tf,
              const typename Block::TensorField<U, R>::IndexRangeType &r,
              typename Block::TensorField<U, R>::DataType *buf)
{
    using IRange = typename Block::TensorField<U, R>::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    // get index space intersection.  The IRange r must be relative to the
    // memory region allocated in buf.
    const IRange indices = r.getIntersection(tf[0].getIndexRange());

    // local block offset
    const MIndex base = indices.getBegin();
    const MIndex offset = base - tf[0].getIndexRange().getBegin();

    // iterate over local (sub) index space and copy into buf
    constexpr size_t Nc = Block::TensorField<U, R>::NComponents;
    for (const auto &p : indices) {
        const size_t i = r.getFlatIndexFromGlobal(base + p);
        const size_t j = tf[0].getIndexRange().getFlatIndex(offset + p);
        for (size_t c = 0; c < Nc; ++c) {
            buf[c + Nc * i] = tf[c][j];
        }
    }
}

/**  @} */
NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* FIELDAOS_H_DOY3MTA8 */
