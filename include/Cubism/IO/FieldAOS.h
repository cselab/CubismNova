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
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

template <Block::FieldClass Class>
struct AOSDriver {
    template <typename Field, typename Range, typename Buffer>
    void write(const Field &, const Range &, Buffer *) const
    {
        throw std::runtime_error(
            "AOSDriver::write not implemented for this Block::FieldClass");
    }
};

template <>
struct AOSDriver<Block::FieldClass::Scalar> {
    template <typename Field, typename Range, typename Buffer>
    void write(const Field &f, const Range &r, Buffer *buf) const
    {
        using MIndex = typename Range::MultiIndex;

        // get index space intersection.  The Range r must be relative to the
        // memory region allocated in buf.
        const Range indices = r.getIntersection(f.getIndexRange());

        // local block offset
        const MIndex base = indices.getBegin();
        const MIndex offset = base - f.getIndexRange().getBegin();

        // iterate over local (sub) index space and copy into buf
        for (const auto &p : indices) {
            const size_t i = r.getFlatIndexFromGlobal(base + p);
            buf[i] = f[offset + p];
        }
    }
};

template <>
struct AOSDriver<Block::FieldClass::Tensor> {
    template <typename Field, typename Range, typename Buffer>
    void write(const Field &f, const Range &r, Buffer *buf) const
    {
        using MIndex = typename Range::MultiIndex;

        // get index space intersection.  The Range r must be relative to the
        // memory region allocated in buf.
        const Range indices = r.getIntersection(f[0].getIndexRange());

        // local block offset
        const MIndex base = indices.getBegin();
        const MIndex offset = base - f[0].getIndexRange().getBegin();

        // iterate over local (sub) index space and copy into buf
        constexpr size_t Nc = Field::NComponents;
        for (const auto &p : indices) {
            const size_t i = r.getFlatIndexFromGlobal(base + p);
            const size_t j = f[0].getIndexRange().getFlatIndex(offset + p);
            for (size_t c = 0; c < Nc; ++c) {
                buf[c + Nc * i] = f[c][j];
            }
        }
    }
};

/**
 * @addtogroup IO
 * @{
 *
 * @brief Write field data into AoS buffer
 * @tparam Buffer Data type of AoS buffer
 * @tparam Field Field type
 * @param f Input field
 * @param r Index space for the copy (describes the memory region of buf)
 * @param buf Output buffer
 *
 * @rst
 * Copy the data from a structure of arrays (SoA) field into an array of
 * structures (AoS) buffer for I/O operation.  This is a low-level function
 * which can be used in high-level I/O interfaces.  The interface is defined for
 * ``Block::FieldClass::Scalar`` and ``Block::FieldClass::Tensor`` fields.  The
 * index range r may describe a sub-region of the index range spanned by the
 * field ``f``.  The size of the output buffer ``buf`` is determined by the
 * index range ``r`` and ``Field::NComponents``.
 * @endrst
 */
template <typename Buffer, typename Field>
void FieldWriteAOS(const Field &f,
                   const typename Field::IndexRangeType &r,
                   Buffer *buf)
{
    AOSDriver<Field::Class> driver;
    driver.write(f, r, buf);
}

// TODO: [fabianw@mavt.ethz.ch; 2020-01-24] Read
/**  @} */

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* FIELDAOS_H_DOY3MTA8 */
