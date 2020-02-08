// File       : Index.h
// Created    : Sun Dec 29 2019 04:56:19 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Defines index space
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef INDEX_H_CPAHO0UM
#define INDEX_H_CPAHO0UM

#include "Cubism/Core/Range.h"
#include <cassert>
#include <cstddef>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Core)

#ifdef CUBISM_32BIT_INDEX
using Index = int;
#else
using Index = std::ptrdiff_t;
#endif

/**
 * @brief Converts indices from and to one dimensional and multi-dimensional,
 * respectively
 * @tparam DIM Dimension of the index space
 */
template <size_t DIM>
struct IndexConverter {
    using MultiIndex = typename Core::Range<Index, DIM>::PointType;

    /**
     * @brief Convert a multi-dimensional index to a one-dimensional index
     * @param p Multi-dimensional index (source)
     * @param extent Extent of index space associated to ``p``
     * @return Flattened index
     */
    size_t getFlatIndex(const MultiIndex &p, const MultiIndex &extent) const
    {
        size_t flat = 0;
        for (size_t i = DIM - 1; i > 0; --i) {
            flat = extent[i - 1] * (p[i] + flat);
        }
        return p[0] + flat;
    }

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-17] remove p.  It must be p(0) for
    // local index
    /**
     * @brief Convert a one-dimensional index to a multi-dimensional index
     * @param i One-dimensional index (source)
     * @param p TODO remove
     * @param extent Extent of index space associated to ``i``
     * @return Multi-dimensional index
     */
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

// forward declaration of iterator type
template <size_t DIM>
class MultiIndexIterator;

/**
 * @brief Rectangular index range
 * @tparam DIM Dimension of the index space
 *
 * Defines a simple consecutive index space.
 */
template <size_t DIM>
class IndexRange : public Core::Range<Index, DIM>
{
public:
    using BaseType = Core::Range<Index, DIM>;
    using typename BaseType::DataType;
    using typename BaseType::PointType;
    using MultiIndex = PointType;

protected:
    using BaseType::begin_;
    using BaseType::end_;
    using BaseType::extent_;

public:
    /** @brief Default constructor (NULL range) */
    IndexRange() : BaseType() {} // NULL range

    /**
     * @brief Construct index range
     * @param e End point (*top right*) of index space. Begin is ``0``.
     *
     * @rst
     * Constructs equal extent in all ``DIM`` dimensions.
     * @endrst
     */
    explicit IndexRange(const DataType e) : BaseType(e) {}
    /**
     * @brief Construct index range
     * @param e End point (*top right*) of index space. Begin is ``0``.
     *
     * @rst
     * Constructs an extent specified the ``DIM``-dimensional ``e``.
     * @endrst
     */
    explicit IndexRange(const PointType &e) : BaseType(e) {}
    /**
     * @brief Construct index range
     * @param b Begin point (*lower left*) of index space.
     * @param e End point (*top right*) of index space.
     *
     * @rst
     * Constructs equal extent in all ``DIM`` dimensions.
     * @endrst
     */
    IndexRange(const DataType b, const DataType e) : BaseType(b, e) {}
    /**
     * @brief Construct index range
     * @param b Begin point (*lower left*) of index space.
     * @param e End point (*top right*) of index space.
     *
     * @rst
     * Constructs an extent specified the ``DIM``-dimensional difference of
     * ``e`` and ``b``.
     * @endrst
     */
    IndexRange(const PointType &b, const PointType &e) : BaseType(b, e) {}

    IndexRange(const IndexRange &c) noexcept = default;
    IndexRange(IndexRange &&c) noexcept = default;
    IndexRange &operator=(const IndexRange &c) = default;
    IndexRange &operator=(IndexRange &&c) = default;

    using iterator = MultiIndexIterator<DIM>;
    iterator begin() noexcept { return iterator(*this, 0); }
    iterator begin() const noexcept { return iterator(*this, 0); }
    iterator end() noexcept { return iterator(*this, this->size()); }
    iterator end() const noexcept { return iterator(*this, this->size()); }
    // iterator end() noexcept { return iterator(*this, this->size() - 1); }
    // iterator end() const noexcept { return iterator(*this, this->size() - 1);
    // }

    /**
     * @brief Get intersection subspace
     * @param o Other index range
     * @return New range for intersection
     */
    IndexRange getIntersection(const IndexRange &o) const
    {
        const auto r = BaseType::getIntersection(o);
        return IndexRange(r.getBegin(), r.getEnd());
    }

    /**
     * @brief Check if index is valid local in this range
     * @param p Local multi-dimensional index
     * @return True if ``p`` is a valid index (exclusive; C-style indexing)
     */
    bool isIndex(const MultiIndex &p) const
    {
        return MultiIndex(0) <= p && p < extent_;
    }

    /**
     * @brief Check if index is valid global in this range
     * @param p Global multi-dimensional index
     * @return True if ``p`` is a valid index (exclusive; C-style indexing)
     */
    bool isGlobalIndex(const MultiIndex &p) const
    {
        return begin_ <= p && p < end_;
    }

    /**
     * @brief Size of index space
     * @return Total number of indices in the index range
     */
    size_t size() const
    {
        if (this->getNullSpace().size() == DIM) {
            return 0;
        }
        return static_cast<size_t>(extent_.prod());
    }

    /**
     * @brief Size of index space
     * @param i Dimension
     * @return  Number of indices along dimension ``i``
     */
    size_t sizeDim(const size_t i) const
    {
        assert(i < DIM);
        return static_cast<size_t>(extent_[i]);
    }

    /**
     * @brief Convert a local multi-dimensional index to a local one-dimensional
     * index
     * @param p Local multi-dimensional index
     * @return Local flattened index
     *
     * @rst
     * Computes a *local* flat index from a *local* multi-dimensional index
     * relative to the index space spanned by this range.
     * @endrst
     */
    size_t getFlatIndex(const MultiIndex &p) const
    {
        assert(MultiIndex(0) <= p && p <= extent_); // inclusive for iterators
        return convert_.getFlatIndex(p, extent_);
    }

    /**
     * @brief Convert a global multi-dimensional index to a local
     * one-dimensional index
     * @param p Global multi-dimensional index
     * @return Local flattened index
     *
     * @rst
     * Computes a *local* flat index from a *global* multi-dimensional index
     * relative to the index space spanned by this range.
     * @endrst
     */
    size_t getFlatIndexFromGlobal(const MultiIndex &p) const
    {
        assert(begin_ <= p && p <= end_); // inclusive for iterators
        return convert_.getFlatIndex(p - begin_, extent_);
    }

    /**
     * @brief Convert a local one-dimensional index to a local multi-dimensional
     * index
     * @param i Local one-dimensional index
     * @return Local multi-dimensional index
     *
     * @rst
     * Computes a *local* multi-dimensional index from a *local* one-dimensional
     * index relative to the index space spanned by this range.
     * @endrst
     */
    MultiIndex getMultiIndex(size_t i) const
    {
        assert(i <= this->size()); // inclusive for iterators
        return convert_.getMultiIndex(i, MultiIndex(0), extent_);
    }

private:
    IndexConverter<DIM> convert_;
};

/**
 * @brief Alias for a multi-dimensional index
 * @tparam DIM Dimension of associated index space
 */
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
protected:
    using IndexRangeType = IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;

public:
    MultiIndexIterator(const IndexRangeType &r, const size_t i) : data_(r, i)
    {
        assert(data_.i <= data_.range.size());
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
        return *this;
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
    const IndexRangeType &getIndexRange() const { return data_.range; }

protected:
    MIIForward<DIM> data_;
};

template <size_t DIM>
class EntityIterator : public MultiIndexIterator<DIM>
{
    using BaseType = MultiIndexIterator<DIM>;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using Entity = Cubism::EntityType;

    using BaseType::data_;

public:
    EntityIterator(const Entity t,
                   const size_t d,
                   const IndexRangeType &r,
                   const size_t i)
        : BaseType(r, i), type_(t), dir_(d)
    {
        assert(dir_ < DIM);
    }
    EntityIterator() = delete;
    EntityIterator(const EntityIterator &c) = default;
    ~EntityIterator() = default;

    EntityIterator &operator=(const EntityIterator &c)
    {
        if (this != &c) {
            BaseType::operator=(c);
            type_ = c.type_;
            dir_ = c.dir_;
        }
        return *this;
    }
    EntityIterator &operator=(const MultiIndex &p)
    {
        BaseType::operator=(p);
        return *this;
    }
    EntityIterator &operator++()
    {
        BaseType::operator++();
        return *this;
    }
    EntityIterator operator++(int)
    {
        EntityIterator tmp(*this);
        ++data_.i;
        data_.forwardMultiIndex();
        return tmp;
    }
    Entity getEntity() const { return type_; }
    size_t getDirection() const { return dir_; }

private:
    Entity type_;
    size_t dir_;
};

NAMESPACE_END(Core)
NAMESPACE_END(Cubism)

#endif /* INDEX_H_CPAHO0UM */
