// File       : StructuredUniform.h
// Created    : Sat Jan 04 2020 11:17:40 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Structured uniform mesh
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef STRUCTUREDUNIFORM_H_L4AIO9HT
#define STRUCTUREDUNIFORM_H_L4AIO9HT

#include "Mesh/StructuredBase.h"

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Mesh)

template <typename TReal, size_t DIM>
class StructuredUniform : public StructuredBase<TReal, DIM>
{
public:
    using BaseMesh = StructuredBase<TReal, DIM>;
    using typename BaseMesh::EntityType;
    using typename BaseMesh::IndexRangeType;
    using typename BaseMesh::MeshClass;
    using typename BaseMesh::MeshHull;
    using typename BaseMesh::MultiIndex;
    using typename BaseMesh::PointType;
    using typename BaseMesh::RangeType;
    using typename BaseMesh::RealType;

private:
    using BaseMesh::crange_;
    using BaseMesh::frange_;
    using BaseMesh::nrange_;
    using BaseMesh::range_;

public:
    StructuredUniform(const PointType &start,
                      const PointType &end,
                      const MultiIndex &cells,
                      const MeshHull type)
        : BaseMesh(start, end, cells, type, MeshClass::Uniform),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    StructuredUniform(const PointType &end,
                      const MultiIndex &cells,
                      const MeshHull type)
        : BaseMesh(end, cells, type, MeshClass::Uniform),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    StructuredUniform(const PointType &gorigin,
                      const RangeType &range,
                      const IndexRangeType &crange,
                      const MeshHull type)
        : BaseMesh(gorigin, range, crange, type, MeshClass::Uniform),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    StructuredUniform(const PointType &gorigin,
                      const RangeType &range,
                      const IndexRangeType &crange,
                      const IndexRangeType &nrange,
                      const std::vector<IndexRangeType> &frange,
                      const MeshHull type)
        : BaseMesh(
              gorigin, range, crange, nrange, frange, type, MeshClass::Uniform),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    StructuredUniform() = delete;
    StructuredUniform(const StructuredUniform &c) = default;
    StructuredUniform(StructuredUniform &&c) noexcept = default;
    StructuredUniform &operator=(const StructuredUniform &c) = delete;
    StructuredUniform &operator=(StructuredUniform &&c) = delete;
    ~StructuredUniform() = default;

protected:
    PointType getCoords_(const MultiIndex &p,
                         const EntityType t,
                         const size_t dir) const override
    {
        if (t == EntityType::Cell) {
            PointType c = PointType(p) + 0.5;
            return range_.getBegin() + c * mesh_spacing_;
        } else if (t == EntityType::Node) {
            PointType c(p);
            return range_.getBegin() + c * mesh_spacing_;
        } else if (t == EntityType::Face) {
            PointType c =
                PointType(p) + 0.5 * (1.0 - PointType::getUnitVector(dir));
            return range_.getBegin() + c * mesh_spacing_;
        } else {
            throw std::runtime_error(
                "StructuredBase::getIndexRange: Unknown entity type t");
        }
        return p;
    }

    RealType getCellVolume_(const MultiIndex &) const override
    {
        return cell_volume_;
    }

    PointType getCellSize_(const MultiIndex &) const override
    {
        return mesh_spacing_;
    }

    PointType getSurface_(const MultiIndex &fi,
                          const MultiIndex &ci,
                          const size_t dir) const override
    {
#ifndef NDEBUG
        assert(dir < DIM);
        for (size_t k = 0; k < fi.size(); ++k) {
            if (k == dir) {
                continue;
            }
            assert(fi[k] - ci[k] == 0);
        }
#endif
        PointType S =
            cell_volume_ / mesh_spacing_[dir] * PointType::getUnitVector(dir);
        S[dir] = (fi[dir] - ci[dir] == 0) ? -S[dir] : S[dir];
        return S;
    }

private:
    const PointType mesh_spacing_;
    const RealType cell_volume_;
};

NAMESPACE_END(Mesh)
NAMESPACE_END(Cubism)

#endif /* STRUCTUREDUNIFORM_H_L4AIO9HT */
