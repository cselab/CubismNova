// File       : Base.h
// Created    : Wed Feb 12 2020 08:23:13 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Base interface for boundary conditions
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef BASE_H_EJIAASP9
#define BASE_H_EJIAASP9

#include "Cubism/Common.h"
#include <cassert>
#include <string>

NAMESPACE_BEGIN(Cubism)
/**
 * @addtogroup BC
 * @{ */
/** @brief Namespace for block field boundary conditions */
NAMESPACE_BEGIN(BC)

/** @brief Boundary information meta data */
struct BoundaryInfo {
    size_t dir;
    size_t side;
    bool is_periodic;
};

/**
 * @brief Boundary condition base class
 * @tparam Lab Type of ``DataLab``
 *
 * @rst
 * Each boundary condition is applied for a specific ``dir < CUBISM_DIMENSION``
 * and corresponding ``side``.
 * @endrst
 * */
template <typename Lab>
class Base
{
    using StencilType = typename Lab::StencilType;

public:
    Base(const size_t dir, const size_t side)
    {
        binfo_.dir = dir;
        binfo_.side = side;
        binfo_.is_periodic = true;
    }
    virtual ~Base() {}

    /**
     * @brief Get boundary information
     * @return ``BoundaryInfo`` structure
     */
    const BoundaryInfo &getBoundaryInfo() const { return binfo_; }

    /**
     * @brief Apply boundary condition
     * @param lab ``DataLab`` where boundary is applied
     */
    virtual void operator()(Lab &) {}

    /**
     * @brief Name of boundary condition
     * @return Name string
     */
    virtual std::string name() const { return std::string("Base"); }

protected:
    BoundaryInfo binfo_;

    bool isValidStencil_(const StencilType &s) const
    {
        assert(binfo_.dir < Lab::IndexRangeType::Dim);
        return !((0 == s.getBegin()[binfo_.dir]) ||
                 (1 == s.getEnd()[binfo_.dir]));
    }
};

NAMESPACE_END(BC)
/**  @} */
NAMESPACE_END(Cubism)

#endif /* BASE_H_EJIAASP9 */
