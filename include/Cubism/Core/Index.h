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

template <size_t DIM>
struct IndexConverter {
    using MultiIndex = typename Core::Range<Index, DIM>::PointType;

    size_t getFlatIndex(const MultiIndex &p, const MultiIndex &extent) const
    {
        size_t flat = 0;
        for (size_t i = DIM - 1; i > 0; --i) {
            flat = extent[i - 1] * (p[i] + flat);
        }
        return p[0] + flat;
    }

    MultiIndex
    getMultiIndex(size_t i, MultiIndex p, const MultiIndex &extent) const
    {
        for (size_t k = 0; k < DIM; ++k) {
            p[k] += i % extent[k];
            i /= extent[k];
        }
        return p;
    }
};

template <>
struct IndexConverter<1> {
    using MultiIndex = typename Core::Range<Index, 1>::PointType;

    size_t getFlatIndex(const MultiIndex &p, const MultiIndex &) const
    {
        return p[0];
    }

    MultiIndex getMultiIndex(size_t i, MultiIndex p, const MultiIndex &) const
    {
        return p += i;
    }
};

template <>
struct IndexConverter<2> {
    using MultiIndex = typename Core::Range<Index, 2>::PointType;

    size_t getFlatIndex(const MultiIndex &p, const MultiIndex &extent) const
    {
        return p[0] + extent[0] * p[1];
    }

    MultiIndex
    getMultiIndex(size_t i, MultiIndex p, const MultiIndex &extent) const
    {
        p[0] += i % extent[0];
        p[1] += i / extent[0];
        return p;
    }
};

template <>
struct IndexConverter<3> {
    using MultiIndex = typename Core::Range<Index, 3>::PointType;

    size_t getFlatIndex(const MultiIndex &p, const MultiIndex &extent) const
    {
        return p[0] + extent[0] * (p[1] + extent[1] * p[2]);
    }

    MultiIndex
    getMultiIndex(size_t i, MultiIndex p, const MultiIndex &extent) const
    {
        p[0] += i % extent[0];
        p[1] += (i / extent[0]) % extent[1];
        p[2] += i / (extent[0] * extent[1]);
        return p;
    }
};

/// @brief Extension of generic Range to index space
template <size_t DIM>
class IndexRange : public Core::Range<Index, DIM>
{
public:
    using BaseType = Core::Range<Index, DIM>;
    using typename BaseType::DataType;
    using typename BaseType::PointType;
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

    /// @brief Return local flattened local MultiIndex
    size_t getFlatIndex(const MultiIndex &p) const
    {
        assert(MultiIndex(0) <= p && p < extent_);
        return convert_.getFlatIndex(p, extent_);
    }

    /// @brief Return local MultiIndex from local flat index
    MultiIndex getMultiIndex(size_t i) const
    {
        assert(i <= this->size());
        return convert_.getMultiIndex(i, MultiIndex(0), extent_);
    }

private:
    PointType extent_;
    IndexConverter<DIM> convert_;
};

template <size_t DIM>
using MultiIndex = typename IndexRange<DIM>::PointType;

template <size_t DIM>
struct MIIForward {
    using IndexRangeType = IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;

    const IndexRangeType &range;
    MultiIndex p;
    size_t i;
    MIIForward(const IndexRangeType &r_, const size_t i_)
        : range(r_), p(r_.getMultiIndex(i_)), i(i_)
    {
    }

    void forwardMultiIndex() { p = range.getMultiIndex(i); }
};

template <>
struct MIIForward<1> {
    using IndexRangeType = IndexRange<1>;
    using MultiIndex = typename IndexRangeType::MultiIndex;

    const IndexRangeType &range;
    MultiIndex p;
    size_t i;
    MIIForward(const IndexRangeType &r_, const size_t i_)
        : range(r_), p(r_.getMultiIndex(i_)), i(i_)
    {
    }

    void forwardMultiIndex() { p[0] = i; }
};

template <>
struct MIIForward<2> {
    using IndexRangeType = IndexRange<2>;
    using MultiIndex = typename IndexRangeType::MultiIndex;

    const IndexRangeType &range;
    const MultiIndex bound;
    MultiIndex p;
    size_t i;
    MIIForward(const IndexRangeType &r_, const size_t i_)
        : range(r_), bound(r_.getExtent() - 1), p(r_.getMultiIndex(i_)), i(i_)
    {
    }

    void forwardMultiIndex()
    {
        p[0] = (p[0] < bound[0]) ? p[0] + 1 : shift_1();
    }
    size_t shift_1()
    {
        p[1] = (p[1] < bound[1]) ? p[1] + 1 : 0;
        return 0;
    }
};

template <>
struct MIIForward<3> {
    using IndexRangeType = IndexRange<3>;
    using MultiIndex = typename IndexRangeType::MultiIndex;

    const IndexRangeType &range;
    const MultiIndex bound;
    MultiIndex p;
    size_t i;
    MIIForward(const IndexRangeType &r_, const size_t i_)
        : range(r_), bound(r_.getExtent() - 1), p(r_.getMultiIndex(i_)), i(i_)
    {
    }

    void forwardMultiIndex()
    {
        p[0] = (p[0] < bound[0]) ? p[0] + 1 : shift_1();
    }
    size_t shift_1()
    {
        p[1] = (p[1] < bound[1]) ? p[1] + 1 : shift_2();
        return 0;
    }
    size_t shift_2()
    {
        p[2] = (p[2] < bound[2]) ? p[2] + 1 : 0;
        return 0;
    }
};

template <size_t DIM>
class MultiIndexIterator
{
    using IndexRangeType = IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using Entity = Cubism::EntityType;

public:
    MultiIndexIterator(const Entity t,
                       const size_t d,
                       const IndexRangeType &r,
                       const size_t i)
        : type_(t), dir_(d), data_(r, i)
    {
        assert(data_.i <= data_.range.size());
        assert(dir_ < DIM);
    }
    MultiIndexIterator() = delete;
    MultiIndexIterator(const MultiIndexIterator &c) = default;
    ~MultiIndexIterator() = default;

    MultiIndexIterator &operator=(const MultiIndexIterator &c)
    {
        if (this != &c) {
            data_.p = c.data_.p;
            data_.i = c.data_.i;
        }
        return *this;
    }
    MultiIndexIterator &operator=(const MultiIndex &p)
    {
        assert(data_.range.contains(p) && (p < data_.range.getEnd()));
        data_.p = p;
        data_.i = data_.range.getFlatIndex(p);
    }
    bool operator==(const MultiIndexIterator &rhs) const
    {
        // relaxed comparison - type_ and dir_ is not considered only index
        return data_.i == rhs.data_.i;
    }
    bool operator!=(const MultiIndexIterator &rhs) const
    {
        // relaxed comparison - type_ and dir_ is not considered only index
        return data_.i != rhs.data_.i;
    }
    MultiIndexIterator &operator++()
    {
        ++data_.i;
        data_.forwardMultiIndex();
        return *this;
    }
    MultiIndexIterator operator++(int)
    {
        MultiIndexIterator tmp(*this);
        ++data_.i;
        data_.forwardMultiIndex();
        return tmp;
    }
    MultiIndex operator-(const MultiIndexIterator &rhs)
    {
        return rhs.data_.p - data_.p;
    }
    const MultiIndex &operator*() const { return data_.p; }
    size_t getFlatIndex() const { return data_.i; }
    MultiIndex getMultiIndex() const { return data_.p; }
    Entity getEntity() const { return type_; }
    size_t getDirection() const { return dir_; }
    const IndexRangeType &getIndexRange() const { return data_.range; }

private:
    const Entity type_;
    const size_t dir_;
    MIIForward<DIM> data_;
};

NAMESPACE_END(Core)
NAMESPACE_END(Cubism)

#endif /* INDEX_H_CPAHO0UM */
