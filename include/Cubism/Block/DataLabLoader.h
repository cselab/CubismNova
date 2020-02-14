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
#include <functional>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

template <typename FieldType, size_t DIM>
struct DataLabLoader {
    using DataType = typename FieldType::DataType;
    using StencilType = Core::Stencil<DIM>;
    using IndexRangeType = typename Core::IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using ID2Field = std::function<const FieldType &(const MultiIndex &)>;
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

// TODO: [fabianw@mavt.ethz.ch; 2020-02-12] specialized versions (notes)
// DataType *dst = block_data_;
// const DataType *src = f0.getData();
// const size_t l0 = range_.sizeDim(0);        // lab stride
// const size_t n0 = curr_range_.sizeDim(0); // data stride
// const size_t bytes0 = n0 * sizeof(DataType);
// const size_t N = curr_range_.size();

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* DATALABLOADER_H_IQCYDR89 */
