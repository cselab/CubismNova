// File       : DataLabLoader.h
// Created    : Fri Feb 14 2020 03:18:12 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Data laboratory load driver specializations
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef DATALABLOADER_H_IQCYDR89
#define DATALABLOADER_H_IQCYDR89

#include "Cubism/Common.h"
#include "Cubism/Core/Index.h"
#include "Cubism/Core/Stencil.h"
#include "Cubism/Core/Vector.h"
#include "Cubism/Math.h"
#include <cstring>
#include <functional>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

template <typename FieldType, size_t DIM>
struct DataLabLoader {
    using DataType = typename FieldType::DataType;
    using StencilType = Core::Stencil<DIM>;
    using IndexRangeType = typename Core::IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using ID2Field = std::function<FieldType &(const MultiIndex &)>;
    using BoolVec = Core::Vector<bool, DIM>;

    StencilType curr_stencil;
    IndexRangeType curr_range;

    void loadInner(const FieldType &src,
                   DataType *dst,
                   const IndexRangeType &rmemory,
                   const MultiIndex &offset)
    {
        const auto range_end = curr_range.end();
        for (auto it = curr_range.begin(); it != range_end; ++it) {
            *(dst + rmemory.getFlatIndex(*it + offset)) =
                src[it.getFlatIndex()];
        }
    }

    void loadGhosts(const MultiIndex &i0,
                    const ID2Field &i2f,
                    DataType *dst,
                    const IndexRangeType &rmemory,
                    const MultiIndex &offset,
                    const BoolVec &periodic,
                    const MultiIndex &skip)
    {
        const MultiIndex one(1);
        const IndexRangeType nbr_range(0, 3);
        const size_t neighbors = nbr_range.size();
        const size_t me = neighbors / 2;
        const MultiIndex curr_extent = curr_range.getExtent();
        const MultiIndex halo_extent =
            curr_extent + curr_stencil.getEnd() - one;
        const MultiIndex stencil_begin = curr_stencil.getBegin();
        for (size_t i = 0; i < neighbors; ++i) {
            if (i == me) {
                continue;
            }
            const MultiIndex bi = nbr_range.getMultiIndex(i) - one;

            typename MultiIndex::DataType isum = 0;
            bool skip_current = false;
            for (size_t j = 0; j < IndexRangeType::Dim; ++j) {
                if (!periodic[j] && bi[j] == skip[j]) {
                    skip_current = true;
                    break;
                }
                isum += Cubism::myAbs(bi[j]);
            }
            if (skip_current) {
                continue;
            }
            if (!curr_stencil.isTensorial() && isum > 1) {
                continue;
            }

            MultiIndex begin;
            MultiIndex end;
            for (size_t j = 0; j < IndexRangeType::Dim; ++j) {
                begin[j] = (bi[j] < 1) ? ((bi[j] < 0) ? stencil_begin[j] : 0)
                                       : curr_extent[j];
                end[j] = (bi[j] < 1) ? ((bi[j] < 0) ? 0 : curr_extent[j])
                                     : halo_extent[j];
            }
            const IndexRangeType halo_range(begin, end);
            const MultiIndex lab_begin = begin + offset;
            const MultiIndex nbr_begin = begin - bi * curr_extent;

            const auto &f = i2f(i0 + bi);
            for (auto &p : halo_range) {
                *(dst + rmemory.getFlatIndex(p + lab_begin)) = f[p + nbr_begin];
            }
        }
    }
};

// specialization for 3D
template <typename FieldType>
struct DataLabLoader<FieldType, 3> {
    using DataType = typename FieldType::DataType;
    using StencilType = Core::Stencil<3>;
    using IndexRangeType = typename Core::IndexRange<3>;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using Index = typename MultiIndex::DataType;
    using ID2Field = std::function<FieldType &(const MultiIndex &)>;
    using BoolVec = Core::Vector<bool, 3>;

    StencilType curr_stencil;
    IndexRangeType curr_range;

    void loadInner(const FieldType &src,
                   DataType *dst,
                   const IndexRangeType &rmemory,
                   const MultiIndex &offset)
    {
        const MultiIndex sstart = curr_stencil.getBegin();
        const MultiIndex extent = curr_range.getExtent();

        const Index sx = -sstart[0];
        const Index sy = -sstart[1];
        const Index ey = sy + (extent[1] / 4) * 4;
        const Index ry = sy + extent[1] - ey;
        const Index sz = -sstart[2];
        const Index ez = sz + extent[2];

        const MultiIndex lextent = rmemory.getExtent(); // lab extent
        const Index lstridex = lextent[0];              // lab stride x
        const Index lslicexy = lextent[0] * lextent[1]; // lab slice xy

        DataType *pdst = dst + rmemory.getFlatIndex(offset + sstart);
        const DataType *psrc = src.getData();
        const Index bstridex = extent[0]; // block stride x
        const size_t bytesx = sizeof(DataType) * bstridex;
        for (Index iz = sz; iz < ez; ++iz) {
            const Index szx = iz * lslicexy + sx;
            for (Index iy = sy; iy < ey; iy += 4) {
                DataType *dst0 = pdst + szx + (iy + 0) * lstridex;
                DataType *dst1 = pdst + szx + (iy + 1) * lstridex;
                DataType *dst2 = pdst + szx + (iy + 2) * lstridex;
                DataType *dst3 = pdst + szx + (iy + 3) * lstridex;
                std::memcpy(dst0, psrc + 0 * bstridex, bytesx);
                std::memcpy(dst1, psrc + 1 * bstridex, bytesx);
                std::memcpy(dst2, psrc + 2 * bstridex, bytesx);
                std::memcpy(dst3, psrc + 3 * bstridex, bytesx);
                psrc += 4 * bstridex;
            }
            if (ry > 0) {
                for (Index iy = ey; iy < ey + ry; ++iy) {
                    DataType *dst0 = pdst + szx + iy * lstridex;
                    std::memcpy(dst0, psrc, bytesx);
                    psrc += bstridex;
                }
            }
        }
    }

    void loadGhosts(const MultiIndex &i0,
                    const ID2Field &i2f,
                    DataType *dst,
                    const IndexRangeType &rmemory,
                    const MultiIndex &offset,
                    const BoolVec &periodic,
                    const MultiIndex &skip)
    {
        const MultiIndex one(1);
        const IndexRangeType nbr_range(0, 3);
        const size_t neighbors = nbr_range.size();
        const size_t me = neighbors / 2;
        const MultiIndex extent = curr_range.getExtent();
        const MultiIndex halo_extent = extent + curr_stencil.getEnd() - one;
        const MultiIndex stencil_begin = curr_stencil.getBegin();
        DataType *pdst = dst + rmemory.getFlatIndex(offset + stencil_begin);
        for (size_t i = 0; i < neighbors; ++i) {
            if (i == me) {
                continue;
            }
            const MultiIndex bi{i % 3 - 1, (i / 3) % 3 - 1, (i / 9) - 1};

            if ((!periodic[0] && bi[0] == skip[0]) ||
                (!periodic[1] && bi[1] == skip[1]) ||
                (!periodic[2] && bi[2] == skip[2])) {
                continue;
            }
            if (!curr_stencil.isTensorial() &&
                (myAbs(bi[0]) + myAbs(bi[1]) + myAbs(bi[2]) > 1)) {
                continue;
            }

            const MultiIndex begin{
                bi[0] < 1 ? (bi[0] < 0 ? stencil_begin[0] : 0) : extent[0],
                bi[1] < 1 ? (bi[1] < 0 ? stencil_begin[1] : 0) : extent[1],
                bi[2] < 1 ? (bi[2] < 0 ? stencil_begin[2] : 0) : extent[2]};

            const MultiIndex end{
                bi[0] < 1 ? (bi[0] < 0 ? 0 : extent[0]) : halo_extent[0],
                bi[1] < 1 ? (bi[1] < 0 ? 0 : extent[1]) : halo_extent[1],
                bi[2] < 1 ? (bi[2] < 0 ? 0 : extent[2]) : halo_extent[2]};

            const auto &f = i2f(i0 + bi);
            const MultiIndex lextent = rmemory.getExtent(); // lab extent
            const Index lstridex = lextent[0];              // lab stride x
            const Index lslicexy = lextent[0] * lextent[1]; // lab slice xy
            const Index sx = begin[0] - stencil_begin[0];
            const size_t bytesx = sizeof(DataType) * (end[0] - begin[0]);
            if (0 == bytesx) {
                continue;
            }
            for (Index iz = begin[2]; iz < end[2]; ++iz) {
                const Index szx = (iz - stencil_begin[2]) * lslicexy + sx;
                if ((end[1] - begin[1]) % 4 != 0) {
                    for (Index iy = begin[1]; iy < end[1]; ++iy) {
                        DataType *dst0 =
                            pdst + szx + (iy - stencil_begin[1]) * lstridex;
                        const DataType *src0 = &f(begin[0] - bi[0] * extent[0],
                                                  iy - bi[1] * extent[1],
                                                  iz - bi[2] * extent[2]);
                        std::memcpy(dst0, src0, bytesx);
                    }
                } else {
                    for (Index iy = begin[1]; iy < end[1]; iy += 4) {
                        DataType *dst0 =
                            pdst + szx + (iy + 0 - stencil_begin[1]) * lstridex;
                        DataType *dst1 =
                            pdst + szx + (iy + 1 - stencil_begin[1]) * lstridex;
                        DataType *dst2 =
                            pdst + szx + (iy + 2 - stencil_begin[1]) * lstridex;
                        DataType *dst3 =
                            pdst + szx + (iy + 3 - stencil_begin[1]) * lstridex;
                        const DataType *src0 = &f(begin[0] - bi[0] * extent[0],
                                                  iy + 0 - bi[1] * extent[1],
                                                  iz - bi[2] * extent[2]);
                        const DataType *src1 = &f(begin[0] - bi[0] * extent[0],
                                                  iy + 1 - bi[1] * extent[1],
                                                  iz - bi[2] * extent[2]);
                        const DataType *src2 = &f(begin[0] - bi[0] * extent[0],
                                                  iy + 2 - bi[1] * extent[1],
                                                  iz - bi[2] * extent[2]);
                        const DataType *src3 = &f(begin[0] - bi[0] * extent[0],
                                                  iy + 3 - bi[1] * extent[1],
                                                  iz - bi[2] * extent[2]);
                        std::memcpy(dst0, src0, bytesx);
                        std::memcpy(dst1, src1, bytesx);
                        std::memcpy(dst2, src2, bytesx);
                        std::memcpy(dst3, src3, bytesx);
                    }
                }
            }
        }
    }
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* DATALABLOADER_H_IQCYDR89 */
