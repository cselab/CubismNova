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

#ifdef CUBISM_32BIT_INDEX
using Index = int;
#else
using Index = std::ptrdiff_t;
#endif

/// @brief Extension of generic Range to index space
template <size_t DIM>
class IndexRange : public Range<Index, DIM>
{
public:
    using BaseType = Range<Index, DIM>;
    using DataType = typename BaseType::DataType;
    using PointType = typename BaseType::PointType;

    // Construction
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

    IndexRange() = delete;
    IndexRange(const IndexRange &c) = default;
    IndexRange(IndexRange &&c) noexcept = default;
    IndexRange &operator=(const IndexRange &c) = default;
    IndexRange &operator=(IndexRange &&c) = default;

    // interface
    size_t size() const { return static_cast<size_t>(extent_.prod()); }

    size_t sizeDim(const size_t i) const
    {
        assert(i < DIM);
        return static_cast<size_t>(extent_[i]);
    }

private:
    PointType extent_;
};

/// @brief Multi-Index type
template <size_t DIM>
using MultiIndex = typename IndexRange<DIM>::PointType;

NAMESPACE_END(Cubism)

#endif /* INDEX_H_CPAHO0UM */
