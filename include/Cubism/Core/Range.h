// File       : Range.h
// Created    : Sun Dec 29 2019 02:03:45 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Generic range class
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef RANGE_H_JDWNYXZA
#define RANGE_H_JDWNYXZA

#include "Core/Vector.h"
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Core)

/**
 * @brief Rectangular range
 * @tparam T Data type that describes coordinates in the range
 * @tparam DIM Dimension
 */
template <typename T, size_t DIM>
class Range
{
public:
    using DataType = T;
    using PointType = Vector<T, DIM>;

    static constexpr size_t Dim = DIM;

    /** @brief Default constructor (NULL range) */
    Range() : begin_(0), end_(0), extent_(end_ - begin_) {}

    /**
     * @brief Construct range
     * @param e End point (*top right*) of rectangle. Begin is ``0``.
     *
     * @rst
     * Constructs equal extent in all ``DIM`` dimensions.
     * @endrst
     */
    explicit Range(const DataType e)
        : begin_(0), end_(e), extent_(end_ - begin_)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }
    /**
     * @brief Construct range
     * @param e End point (*top right*) of rectangle. Begin is ``0``.
     *
     * @rst
     * Constructs an extent specified the ``DIM``-dimensional ``e``.
     * @endrst
     */
    explicit Range(const PointType &e)
        : begin_(0), end_(e), extent_(end_ - begin_)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }
    /**
     * @brief Construct range
     * @param b Begin point (*lower left*) of rectangle.
     * @param e End point (*top right*) of rectangle.
     *
     * @rst
     * Constructs equal extent in all ``DIM`` dimensions.
     * @endrst
     */
    Range(const DataType b, const DataType e)
        : begin_(b), end_(e), extent_(end_ - begin_)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }
    /**
     * @brief Construct range
     * @param b Begin point (*lower left*) of rectangle.
     * @param e End point (*top right*) of rectangle.
     *
     * @rst
     * Constructs an extent specified the ``DIM``-dimensional difference of
     * ``e`` and ``b``.
     * @endrst
     */
    Range(const PointType &b, const PointType &e)
        : begin_(b), end_(e), extent_(end_ - begin_)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }

    Range(const Range &c) = default;
    Range(Range &&c) noexcept = default;
    Range &operator=(const Range &c) = default;
    Range &operator=(Range &&c) = default;
    virtual ~Range() = default;

    /**
     * @brief Set range begin
     * @param b New begin point
     */
    void setBegin(const PointType &b)
    {
        begin_ = b;
        extent_ = end_ - begin_;
        check_("RangeSetBegin");
    }

    /**
     * @brief Set range end
     * @param e New end point
     */
    void setEnd(const PointType &e)
    {
        end_ = e;
        extent_ = end_ - begin_;
        check_("RangeSetEnd");
    }

    /**
     * @brief Get range begin
     * @return Begin point
     */
    PointType getBegin() const { return begin_; }

    /**
     * @brief Get range end
     * @return End point
     */
    PointType getEnd() const { return end_; }

    /**
     * @brief Get range extent
     * @return Range extent
     */
    PointType getExtent() const { return extent_; }

    /**
     * @brief Get range volume
     * @return Range volume
     */
    DataType getVolume() const { return extent_.prod(); }

    /**
     * @brief Check if other range is contained in this range
     * @param o Other range
     * @return True if ``o`` is contained (inclusive)
     */
    bool isContained(const Range &o) const
    {
        return begin_ <= o.begin_ && o.end_ <= end_;
    }

    /**
     * @brief Check if point is contained in this range
     * @param p Point
     * @return True if ``p`` is contained (inclusive)
     */
    bool isContained(const PointType &p) const
    {
        return begin_ <= p && p <= end_;
    }

    /**
     * @brief Check if other range intersects this range
     * @param o Other range
     * @return True if ``o`` intersects this range
     */
    bool isIntersecting(const Range &o) const
    {
        return begin_ < o.end_ && o.begin_ < end_;
    }

    /**
     * @brief Check if other range is equal to this range
     * @param o Other range
     */
    bool operator==(const Range &o) const
    {
        return begin_ == o.begin_ && end_ == o.end_;
    }
    /**
     * @brief Check if other range is not equal to this range
     * @param o Other range
     */
    bool operator!=(const Range &o) const { return !(*this == o); }

protected:
    PointType begin_, end_, extent_;

    void check_(const std::string &where) const
    {
        if (begin_ > end_) {
            throw std::runtime_error(where +
                                     ": begin_ must be smaller than end_");
        }
    }
};

template <typename T, size_t DIM>
constexpr size_t Range<T, DIM>::Dim;

NAMESPACE_END(Core)
NAMESPACE_END(Cubism)

#endif /* RANGE_H_JDWNYXZA */
