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
    using STDFunction = std::function<FieldType &(const MultiIndex &)>;
    using BoolVec = Core::Vector<bool, DIM>;

    StencilType curr_stencil;
    IndexRangeType curr_range;
    IndexRangeType curr_labrange;

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

    template <typename Functor = STDFunction>
    void loadGhosts(const MultiIndex &i0,
                    const Functor &i2f,
                    DataType *dst,
                    const IndexRangeType &rmemory,
                    const MultiIndex &offset,
                    const BoolVec &periodic,
                    const MultiIndex &skip)
    {
        const IndexRangeType nbr_range(0, 3);
        const size_t neighbors = nbr_range.size();
        const size_t me = neighbors / 2;
        const MultiIndex curr_extent = curr_range.getExtent();
        const MultiIndex halo_extent = curr_extent + curr_stencil.getEnd() - 1;
        const MultiIndex stencil_begin = curr_stencil.getBegin();
        for (size_t i = 0; i < neighbors; ++i) {
            if (i == me) {
                continue;
            }
            const MultiIndex bi = nbr_range.getMultiIndex(i) - 1;

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

            const auto &f = i2f(i0 + bi); // neighbor block field
            const MultiIndex nbr_extent = f.getIndexRange().getExtent();
            MultiIndex begin;
            MultiIndex end;
            MultiIndex shift(0);
            for (size_t j = 0; j < IndexRangeType::Dim; ++j) {
                if (bi[j] < 1) {
                    if (bi[j] < 0) {
                        shift[j] = nbr_extent[j] - curr_extent[j];
                        begin[j] = stencil_begin[j];
                        end[j] = 0;
                    } else {
                        begin[j] = 0;
                        end[j] = curr_extent[j];
                    }
                } else {
                    begin[j] = curr_extent[j];
                    end[j] = halo_extent[j];
                }
            }
            const IndexRangeType halo_range(begin, end);
            const MultiIndex lab_begin = begin + offset;
            const MultiIndex nbr_begin = begin - bi * curr_extent + shift;
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
    using STDFunction = std::function<FieldType &(const MultiIndex &)>;
    using BoolVec = Core::Vector<bool, 3>;

    StencilType curr_stencil;
    IndexRangeType curr_range;
    IndexRangeType curr_labrange;

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

    template <typename Functor = STDFunction>
    void loadGhosts(const MultiIndex &i0,
                    Functor &i2f,
                    DataType *dst,
                    const IndexRangeType &rmemory,
                    const MultiIndex &offset,
                    const BoolVec &periodic,
                    const MultiIndex &skip)
    {
        const IndexRangeType nbr_range(0, 3);
        const size_t neighbors = nbr_range.size();
        const size_t me = neighbors / 2;
        const MultiIndex extent = curr_range.getExtent();
        const MultiIndex halo_extent = extent + curr_stencil.getEnd() - 1;
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

            const auto &f = i2f(i0 + bi);
            const MultiIndex nbr_extent = f.getIndexRange().getExtent();
            const MultiIndex begin{
                bi[0] < 1 ? (bi[0] < 0 ? stencil_begin[0] : 0) : extent[0],
                bi[1] < 1 ? (bi[1] < 0 ? stencil_begin[1] : 0) : extent[1],
                bi[2] < 1 ? (bi[2] < 0 ? stencil_begin[2] : 0) : extent[2]};
            const MultiIndex end{
                bi[0] < 1 ? (bi[0] < 0 ? 0 : extent[0]) : halo_extent[0],
                bi[1] < 1 ? (bi[1] < 0 ? 0 : extent[1]) : halo_extent[1],
                bi[2] < 1 ? (bi[2] < 0 ? 0 : extent[2]) : halo_extent[2]};
            const MultiIndex shift{bi[0] < 0 ? nbr_extent[0] - extent[0] : 0,
                                   bi[1] < 0 ? nbr_extent[1] - extent[1] : 0,
                                   bi[2] < 0 ? nbr_extent[2] - extent[2] : 0};
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
                        const DataType *src0 =
                            &f(begin[0] - bi[0] * extent[0] + shift[0],
                               iy - bi[1] * extent[1] + shift[1],
                               iz - bi[2] * extent[2] + shift[2]);
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
                        const DataType *src0 =
                            &f(begin[0] - bi[0] * extent[0] + shift[0],
                               iy + 0 - bi[1] * extent[1] + shift[1],
                               iz - bi[2] * extent[2] + shift[2]);
                        const DataType *src1 =
                            &f(begin[0] - bi[0] * extent[0] + shift[0],
                               iy + 1 - bi[1] * extent[1] + shift[1],
                               iz - bi[2] * extent[2] + shift[2]);
                        const DataType *src2 =
                            &f(begin[0] - bi[0] * extent[0] + shift[0],
                               iy + 2 - bi[1] * extent[1] + shift[1],
                               iz - bi[2] * extent[2] + shift[2]);
                        const DataType *src3 =
                            &f(begin[0] - bi[0] * extent[0] + shift[0],
                               iy + 3 - bi[1] * extent[1] + shift[1],
                               iz - bi[2] * extent[2] + shift[2]);
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

template <typename FContainer, Cubism::FieldClass Class, size_t RANK>
struct ScalarFieldMap {
    using BaseType = typename FContainer::BaseType;
    FContainer &fields;
    ScalarFieldMap(FContainer &f) : fields(f) {}
    typename BaseType::FieldType &
    operator()(const size_t i, const size_t c, const size_t)
    {
        assert(c < BaseType::NComponents);
        return fields[i][c];
    }
    const typename BaseType::FieldType &
    operator()(const size_t i, const size_t c, const size_t) const
    {
        assert(c < BaseType::NComponents);
        return fields[i][c];
    }
};

template <typename FContainer, Cubism::FieldClass Class>
struct ScalarFieldMap<FContainer, Class, 0> {
    using BaseType = typename FContainer::BaseType;
    FContainer &fields;
    ScalarFieldMap(FContainer &f) : fields(f) {}
    typename BaseType::FieldType &
    operator()(const size_t i, const size_t, const size_t)
    {
        return fields[i];
    }
    const typename BaseType::FieldType &
    operator()(const size_t i, const size_t, const size_t) const
    {
        return fields[i];
    }
};

template <typename FContainer, size_t RANK>
struct ScalarFieldMap<FContainer, Cubism::FieldClass::FaceContainer, RANK> {
    using BaseType = typename FContainer::BaseType;
    FContainer &fields;
    ScalarFieldMap(FContainer &f) : fields(f) {}
    typename BaseType::FieldType &
    operator()(const size_t i, const size_t c, const size_t d)
    {
        assert(c < BaseType::NComponents);
        assert(d < BaseType::IndexRangeType::Dim);
        return fields[i][d][c];
    }
    const typename BaseType::FieldType &
    operator()(const size_t i, const size_t c, const size_t d) const
    {
        assert(c < BaseType::NComponents);
        assert(d < BaseType::IndexRangeType::Dim);
        return fields[i][d][c];
    }
};

template <typename FContainer>
struct ScalarFieldMap<FContainer, Cubism::FieldClass::FaceContainer, 0> {
    using BaseType = typename FContainer::BaseType;
    FContainer &fields;
    ScalarFieldMap(FContainer &f) : fields(f) {}
    typename BaseType::FieldType &
    operator()(const size_t i, const size_t, const size_t d)
    {
        assert(d < BaseType::IndexRangeType::Dim);
        return fields[i][d];
    }
    const typename BaseType::FieldType &
    operator()(const size_t i, const size_t, const size_t d) const
    {
        assert(d < BaseType::IndexRangeType::Dim);
        return fields[i][d];
    }
};

template <typename FContainer, Cubism::FieldClass Class, size_t RANK>
class PeriodicIndexFunctor
{
public:
    using ScalarField = typename FContainer::BaseType::FieldType;
    using IndexRangeType = typename ScalarField::IndexRangeType;
    using MultiIndex = typename IndexRangeType::MultiIndex;

    PeriodicIndexFunctor(FContainer &fields,
                         const IndexRangeType &range,
                         const size_t comp = 0,
                         const size_t fdir = 0)
        : fields_(fields), range_(range), extent_(range_.getExtent()),
          comp_(comp), face_dir_(fdir)
    {
    }
    PeriodicIndexFunctor() = delete;
    PeriodicIndexFunctor(const PeriodicIndexFunctor &c) = default;
    PeriodicIndexFunctor(PeriodicIndexFunctor &&c) = default;
    PeriodicIndexFunctor &operator=(const PeriodicIndexFunctor &c) = default;
    PeriodicIndexFunctor &operator=(PeriodicIndexFunctor &&c) = default;

    ScalarField &operator()(const MultiIndex &p)
    {
        return fields_(range_.getFlatIndex(periodic_(p)), comp_, face_dir_);
    }

    const ScalarField &operator()(const MultiIndex &p) const
    {
        return fields_(range_.getFlatIndex(periodic_(p)), comp_, face_dir_);
    }

private:
    ScalarFieldMap<FContainer, Class, RANK> fields_;
    const IndexRangeType range_;
    const MultiIndex extent_;
    const size_t comp_;     // component
    const size_t face_dir_; // face direction

    MultiIndex periodic_(MultiIndex p) const
    {
        for (size_t i = 0; i < IndexRangeType::Dim; ++i) {
            p[i] = (p[i] + extent_[i]) % extent_[i];
        }
        return p;
    }
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* DATALABLOADER_H_IQCYDR89 */
