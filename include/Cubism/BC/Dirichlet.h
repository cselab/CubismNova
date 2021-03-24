// File       : Dirichlet.h
// Created    : Sat Feb 15 2020 01:52:58 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Dirichlet boundary conditions
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef DIRICHLET_H_OIZZ1SS4
#define DIRICHLET_H_OIZZ1SS4

#include "Cubism/BC/Base.h"
#include <cassert>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(BC)

/**
 * @brief Dirichlet BC
 * @tparam Lab Type of ``FieldLab``
 *
 * @rst
 * Constant value Dirichlet boundary condition
 * @endrst
 * */
template <typename Lab>
class Dirichlet : public BC::Base<Lab>
{
    using BaseType = BC::Base<Lab>;
    using BaseType::binfo_;

    using DataType = typename Lab::DataType;

public:
    /**
     * @brief Main constructor
     * @param dir Direction in which to apply the boundary
     * @param side On which side along direction ``dir``
     * @param val Boundary value
     */
    Dirichlet(const size_t dir, const size_t side, const DataType &val)
        : BaseType(dir, side), value_(val)
    {
        binfo_.is_periodic = false;
    }

    /**
     * @brief Apply boundary condition
     * @param lab Lab on which the boundary is applied
     */
    void operator()(Lab &lab) override { apply_(lab); }

    /**
     * @brief Name of boundary condition
     * @return Name string
     */
    std::string name() const override { return std::string("Dirichlet"); }

    /**
     * @brief Get boundary value
     * @return Reference to ``DataType``
     */
    DataType &getValue() { return value_; }
    /**
     * @brief Get boundary value
     * @return ``const`` reference to ``DataType``
     */
    const DataType &getValue() const { return value_; }

private:
    using IndexRangeType = typename Lab::IndexRangeType;
    using MultiIndex = typename Lab::MultiIndex;
    using Index = typename MultiIndex::DataType;
    using Stencil = typename Lab::StencilType;

    DataType value_;

    void apply_(Lab &lab) const
    {
        assert(binfo_.dir < IndexRangeType::Dim);
        assert(0 == binfo_.side || 1 == binfo_.side);

        const Stencil &stencil = lab.getActiveStencil();
        if (!this->isValidStencil_(stencil)) {
            return; // nothing to do; zero stencil width for binfo_.dir
        }

        MultiIndex extent;
        if (stencil.isTensorial()) {
            extent = lab.getActiveLabRange().getExtent();
        } else {
            extent = lab.getActiveRange().getExtent();
        }

        const MultiIndex sbegin = stencil.getBegin();
        const MultiIndex send = stencil.getEnd();
        MultiIndex start(0);
        if (0 == binfo_.side) {
            extent[binfo_.dir] = -sbegin[binfo_.dir];
            if (stencil.isTensorial()) {
                start = sbegin;
            } else {
                start[binfo_.dir] = sbegin[binfo_.dir];
            }
        } else {
            extent[binfo_.dir] = send[binfo_.dir] - 1;
            if (stencil.isTensorial()) {
                start = sbegin;
            }
            start[binfo_.dir] = lab.getActiveRange().getExtent()[binfo_.dir];
        }
        const IndexRangeType slab(extent);
        for (const auto &p : slab) {
            lab[p + start] = value_;
        }
    }
};

NAMESPACE_END(BC)
NAMESPACE_END(Cubism)

#endif /* DIRICHLET_H_OIZZ1SS4 */
