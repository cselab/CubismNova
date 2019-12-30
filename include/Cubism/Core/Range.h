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

template <typename T, size_t DIM>
class Range
{
public:
    using DataType = T;
    using PointType = Cubism::Vector<T, DIM>;

    static constexpr size_t Dim = DIM;

    // Construction
    explicit Range(const DataType e) : begin_(), end_(e)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }
    explicit Range(const PointType &e) : begin_(), end_(e)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }
    Range(const DataType b, const DataType e) : begin_(b), end_(e)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }
    Range(const PointType &b, const PointType &e) : begin_(b), end_(e)
    {
        static_assert(DIM > 0, "Spatial dimension DIM must be greater than 0");
        check_("RangeConstruction");
    }

    Range() = delete;
    Range(const Range &c) = default;
    Range(Range &&c) noexcept = default;
    Range &operator=(const Range &c) = default;
    Range &operator=(Range &&c) = default;
    virtual ~Range() = default;

    /// @brief Set range begin
    ///
    /// @param b New begin
    void setBegin(const PointType &b)
    {
        begin_ = b;
        check_("RangeSetBegin");
    }

    /// @brief Set range end
    ///
    /// @param e New end
    void setEnd(const PointType &e)
    {
        end_ = e;
        check_("RangeSetEnd");
    }

    /// @brief Get current range begin
    ///
    /// @return PointType begin
    PointType getBegin() const { return begin_; }

    /// @brief Get current range end
    ///
    /// @return PointType end
    PointType getEnd() const { return end_; }

    /// @brief Get current range extent
    ///
    /// @return PointType range extent
    PointType getExtent() const { return end_ - begin_; }

    /// @brief Get range volume
    ///
    /// @return DataType range volume
    DataType getVolume() const { return getExtent().prod(); }

    /// @brief Check if other range is contained in this range
    ///
    /// @param o Other range
    bool contains(const Range &o) const
    {
        return begin_ <= o.begin_ && o.end_ <= end_;
    }

    /// @brief Check if point is contained in this range
    ///
    /// @param p Point
    bool contains(const PointType &p) const { return begin_ <= p && p <= end_; }

    /// @brief Check if other range intersects this range
    ///
    /// @param o Other range
    bool intersect(const Range &o) const
    {
        return begin_ < o.end_ || o.begin_ < end_;
    }

    /// @brief Check if other range is equal to this range
    ///
    /// @param o Other range
    bool operator==(const Range &o) const
    {
        return begin_ == o.begin_ && end_ == o.end_;
    }
    bool operator!=(const Range &o) const { return !(*this == o); }

protected:
    PointType begin_, end_;

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

NAMESPACE_END(Cubism)

#endif /* RANGE_H_JDWNYXZA */
