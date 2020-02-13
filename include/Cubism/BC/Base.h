// File       : Base.h
// Created    : Wed Feb 12 2020 08:23:13 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Base interface for boundary conditions
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef BASE_H_EJIAASP9
#define BASE_H_EJIAASP9

#include "Cubism/Common.h"

NAMESPACE_BEGIN(Cubism)
/**
 * @addtogroup BC
 * @{ */
/** @brief Namespace for block field boundary conditions */
NAMESPACE_BEGIN(BC)

/** @brief Boundary information meta data */
struct BoundaryInfo {
    bool is_periodic;
    size_t dir;
    size_t side;
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
public:
    Base() = default;
    virtual ~Base() {}

    /**
     * @brief Apply boundary condition
     * @param lab ``DataLab`` where boundary is applied
     */
    virtual void operator()(Lab &) {}

    /**
     * @brief Get boundary information
     * @return ``BoundaryInfo`` structure
     */
    BoundaryInfo getBoundaryInfo() const { return binfo_; }

protected:
    BoundaryInfo binfo_;
};

NAMESPACE_END(BC)
/**  @} */
NAMESPACE_END(Cubism)

#endif /* BASE_H_EJIAASP9 */
