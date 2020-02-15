// File       : Absorbing.h
// Created    : Fri Feb 14 2020 12:06:10 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Absorbing boundary conditions
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef ABSORBING_H_YUEH2NIZ
#define ABSORBING_H_YUEH2NIZ

#include "Cubism/BC/Base.h"
#include <cassert>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(BC)

/**
 * @brief Zeroth-Order absorbing BC
 * @tparam Lab Type of ``DataLab``
 *
 * @rst
 * Zeroth-Order absorbing boundary condition
 * @endrst
 * */
template <typename Lab>
class Absorbing : public BC::Base<Lab>
{
    using BaseType = BC::Base<Lab>;
    using BaseType::binfo_;

public:
    /**
     * @brief Main constructor
     * @param dir Direction in which to apply the boundary
     * @param side On which side along direction ``dir``
     */
    Absorbing(const size_t dir, const size_t side) : BaseType(dir, side)
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
    std::string name() const override
    {
        return std::string("Zeroth-Order Absorbing");
    }

private:
    using DataType = typename Lab::DataType;
    using IndexRangeType = typename Lab::IndexRangeType;
    using MultiIndex = typename Lab::MultiIndex;
    using Index = typename MultiIndex::DataType;
    using Stencil = typename Lab::StencilType;

    void apply_(Lab &lab) const
    {
        assert(binfo_.dir < IndexRangeType::Dim);
        assert(0 == binfo_.side || 1 == binfo_.side);

        const Stencil &stencil = lab.getActiveStencil();
        if (!this->isValidStencil_(stencil)) {
            return; // nothing to do; zero stencil width for binfo_.dir
        }
        const MultiIndex extent = lab.getActiveRange().getExtent();
        const MultiIndex src = (0 == binfo_.side)
                                   ? MultiIndex::getUnitVector(binfo_.dir)
                                   : -MultiIndex::getUnitVector(binfo_.dir);
        const Index sdir = (0 == binfo_.side) ? -1 : 1;

        MultiIndex slice;
        MultiIndex begin(0);
        Index end;
        if (stencil.isTensorial()) {
            slice = lab.getActiveLabRange().getExtent();
            begin = stencil.getBegin();
        } else {
            slice = lab.getActiveRange().getExtent();
        }
        if (0 == binfo_.side) {
            begin[binfo_.dir] = -1;
            end = -stencil.getBegin()[binfo_.dir];
        } else {
            begin[binfo_.dir] = extent[binfo_.dir];
            end = stencil.getEnd()[binfo_.dir] - 1;
        }
        slice[binfo_.dir] = 1;
        const IndexRangeType srange(slice);
        for (const auto &p : srange) {
            MultiIndex q = p + begin;
            const DataType val = lab[q + src];
            lab[q] = val;
            for (Index i = 1; i < end; ++i) {
                q[binfo_.dir] += sdir;
                lab[q] = val;
            }
        }
    }
};

NAMESPACE_END(BC)
NAMESPACE_END(Cubism)

#endif /* ABSORBING_H_YUEH2NIZ */
