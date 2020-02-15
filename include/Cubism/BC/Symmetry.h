// File       : Symmetry.h
// Created    : Sat Feb 15 2020 03:22:19 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Symmetry boundary condition
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef SYMMETRY_H_8HROEJBQ
#define SYMMETRY_H_8HROEJBQ

#include "Cubism/BC/Base.h"
#include <cassert>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(BC)

/**
 * @brief Symmetry BC (reflecting)
 * @tparam Lab Type of ``DataLab``
 *
 * @rst
 * Symmetry/reflecting boundary condition
 * @endrst
 * */
template <typename Lab>
class Symmetry : public BC::Base<Lab>
{
    using BaseType = BC::Base<Lab>;
    using BaseType::binfo_;

    using DataType = typename Lab::DataType;

public:
    /**
     * @brief Main constructor
     * @param dir Direction in which to apply the boundary
     * @param side On which side along direction ``dir``
     */
    Symmetry(const size_t dir, const size_t side, const DataType sign = 1)
        : BaseType(dir, side), sign_(sign)
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
    std::string name() const override { return std::string("Symmetry"); }

private:
    using IndexRangeType = typename Lab::IndexRangeType;
    using MultiIndex = typename Lab::MultiIndex;
    using Index = typename MultiIndex::DataType;
    using Stencil = typename Lab::StencilType;

    const DataType sign_;

    void apply_(Lab &lab) const
    {
        assert(binfo_.dir < IndexRangeType::Dim);
        assert(0 == binfo_.side || 1 == binfo_.side);

        const Stencil &stencil = lab.getActiveStencil();
        if (!this->isValidStencil_(stencil)) {
            return; // zero stencil width for binfo_.dir
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
        Index roffset = 0;
        if (0 == binfo_.side) {
            extent[binfo_.dir] = -sbegin[binfo_.dir];
            roffset = -1;
            if (stencil.isTensorial()) {
                start = sbegin;
            } else {
                start[binfo_.dir] = sbegin[binfo_.dir];
            }
        } else {
            extent[binfo_.dir] = send[binfo_.dir] - 1;
            const MultiIndex N = lab.getActiveRange().getExtent();
            roffset = 2 * N[binfo_.dir] - 1;
            if (stencil.isTensorial()) {
                start = sbegin;
            }
            start[binfo_.dir] = N[binfo_.dir];
        }
        const IndexRangeType slab(extent);
        for (const auto &p : slab) {
            const MultiIndex q = p + start;
            MultiIndex r(q);
            r[binfo_.dir] = roffset - r[binfo_.dir];
            lab[q] = sign_ * lab[r];
        }
    }
};

NAMESPACE_END(BC)
NAMESPACE_END(Cubism)

#endif /* SYMMETRY_H_8HROEJBQ */
