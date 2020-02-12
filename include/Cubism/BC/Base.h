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

// TODO: [fabianw@mavt.ethz.ch; 2020-02-12] doxy
template <typename Lab>
class Base
{
public:
    Base() = default;
    virtual ~Base() {}

    // XXX: [fabianw@mavt.ethz.ch; 2020-02-12] pure virtual ?
    virtual void operator()(Lab &, const double = 0) {}
    virtual bool isPeriodic() const { return true; }
};

NAMESPACE_END(BC)
/**  @} */
NAMESPACE_END(Cubism)

#endif /* BASE_H_EJIAASP9 */
