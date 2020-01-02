// File       : Index.h
// Created    : Sun Dec 29 2019 04:56:19 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Defines index space
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef INDEX_H_CPAHO0UM
#define INDEX_H_CPAHO0UM

#include "Core/Range.h"
#include <cassert>
#include <cstddef>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Core)

#ifdef CUBISM_32BIT_INDEX
using Index = int;
#else
using Index = std::ptrdiff_t;
#endif

/// @brief Extension of generic Range to index space
template <size_t DIM>
class IndexRange : public Core::Range<Index, DIM>
{
public:
    using BaseType = Range<Index, DIM>;
    using DataType = typename BaseType::DataType;
    using PointType = typename BaseType::PointType;
    using MultiIndex = PointType;

    // Construction
    IndexRange() : BaseType(0), extent_(0) {} // NULL range

    explicit IndexRange(const DataType e)
        : BaseType(e), extent_(this->end_ - this->begin_)
    {
    }
    explicit IndexRange(const PointType &e)
        : BaseType(e), extent_(this->end_ - this->begin_)
    {
    }
    IndexRange(const DataType b, const DataType e)
        : BaseType(b, e), extent_(this->end_ - this->begin_)
    {
    }
    IndexRange(const PointType &b, const PointType &e)
        : BaseType(b, e), extent_(this->end_ - this->begin_)
    {
    }

    IndexRange(const IndexRange &c) = default;
    IndexRange(IndexRange &&c) noexcept = default;
    IndexRange &operator=(const IndexRange &c) = default;
    IndexRange &operator=(IndexRange &&c) = default;

    /// @brief Get current range extent
    ///
    /// @return PointType range extent
    PointType getExtent() const { return extent_; }

    /// @brief Return size of this space (number of indices)
    size_t size() const { return static_cast<size_t>(extent_.prod()); }

    /// @brief Return number of indices along dimension i
    size_t sizeDim(const size_t i) const
    {
        assert(i < DIM);
        return static_cast<size_t>(extent_[i]);
    }

    /// @brief Return flattened MultiIndex
    size_t getFlatIndex(const MultiIndex &p) const
    {
        assert(MultiIndex(0) <= p && p < extent_);
        size_t flat = 0;
        for (size_t i = DIM - 1; i > 0; --i) {
            flat = extent_[i - 1] * (p[i] + flat);
        }
        return p[0] + flat;
    }

    /// @brief Return MultiIndex from flat index
    MultiIndex getMultiIndex(size_t i) const
    {
        assert(i < this->size());
        MultiIndex p = this->begin_;
        for (size_t k = 0; k < DIM; ++k) {
            p[k] += i % extent_[k];
            i /= extent_[k];
        }
        return p;
    }

private:
    PointType extent_;
};

NAMESPACE_END(Core)
NAMESPACE_END(Cubism)

#endif /* INDEX_H_CPAHO0UM */
