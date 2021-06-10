// File       : Stencil.h
// Created    : Tue Feb 11 2020 04:55:37 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Type for stencil description
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef STENCIL_H_WOXG1TZC
#define STENCIL_H_WOXG1TZC

#include "Cubism/Common.h"
#include "Cubism/Core/Index.h"
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Core)

/**
 * @brief Describes a stencil
 * @tparam DIM Stencil dimensionality
 * */
template <size_t DIM = CUBISM_DIMENSION>
class Stencil
{
    using BaseType = Core::IndexRange<DIM>;

public:
    using DataType = typename BaseType::DataType;
    using PointType = typename BaseType::PointType;
    using MultiIndex = PointType;

    /** @brief Default constructor */
    Stencil() : begin_(0), end_(1), is_tensorial_(false) {}

    /**
     * @brief Main constructor
     * @param b Begin of stencil
     * @param e End of stencil
     * @param tensorial Flag for tensorial stencil type
     *
     * @rst
     * A symmetric stencil from -1 to +1 is constructed by ``Stencil(-1,2)``.
     * The stencil end is exclusive.
     * @endrst
     */
    Stencil(const DataType b, const DataType e, const bool tensorial = false)
        : begin_(b), end_(e), is_tensorial_(tensorial)
    {
        check_();
    }

    /**
     * @brief Main constructor
     * @param b Begin of stencil
     * @param e End of stencil
     * @param tensorial Flag for tensorial stencil type
     *
     * @rst
     * An arbitrary stencil in 2D may be constructed as follows
     * ``Stencil<2>({-1, -2}, {1, 3})``. The stencil end is exclusive.
     * @endrst
     */
    Stencil(const PointType &b,
            const PointType &e,
            const bool tensorial = false)
        : begin_(b), end_(e), is_tensorial_(tensorial)
    {
        check_();
    }

    /**
     * @brief Get stencil begin
     * @return Multi-dimensional index of stencil begin
     *
     * The stencil begin is inclusive.
     */
    MultiIndex getBegin() const { return begin_; }

    /**
     * @brief Get stencil end
     * @return Multi-dimensional index of stencil end
     *
     * The stencil end is exclusive.
     */
    MultiIndex getEnd() const { return end_; }

    /**
     * @brief Check if stencil is tensorial
     * @return True if tensorial
     *
     * Tensorial stencils take into account edge and corner ghosts.  A
     * non-tensorial stencil accounts for face ghosts only.
     */
    bool isTensorial() const { return is_tensorial_; }

private:
    MultiIndex begin_; // <= 0 (inclusive)
    MultiIndex end_;   // > 0 (exclusive)
    bool is_tensorial_;

    void check_() const
    {
        if (begin_ > MultiIndex(0)) {
            throw std::runtime_error("Stencil: begin_ must be <= 0");
        }
        if (end_ < MultiIndex(1)) {
            throw std::runtime_error("Stencil: end_ must be > 0");
        }
    }
};

NAMESPACE_END(Core)
NAMESPACE_END(Cubism)

#endif /* STENCIL_H_WOXG1TZC */
