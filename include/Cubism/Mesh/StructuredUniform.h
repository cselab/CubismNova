// File       : StructuredUniform.h
// Created    : Sat Jan 04 2020 11:17:40 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Structured uniform mesh
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef STRUCTUREDUNIFORM_H_L4AIO9HT
#define STRUCTUREDUNIFORM_H_L4AIO9HT

#include "Cubism/Math.h"
#include "Cubism/Mesh/StructuredBase.h"
#include <limits>
#include <memory>
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Mesh)
/**
 * @defgroup Mesh Mesh Classes
 * The members of this group belong to mesh classes.
 */

/**
 * @ingroup Mesh
 * @brief Structured uniform mesh
 * @tparam TReal Float type for mesh entities
 * @tparam DIM Mesh dimension
 * */
template <typename TReal, size_t DIM>
class StructuredUniform : public StructuredBase<TReal, DIM>
{
public:
    /** @brief Base mesh type */
    using BaseMesh = StructuredBase<TReal, DIM>;
    using typename BaseMesh::EntityType;
    using typename BaseMesh::IndexRangeType;
    using typename BaseMesh::MeshIntegrity;
    using typename BaseMesh::MultiIndex;
    using typename BaseMesh::PointType;
    using typename BaseMesh::RangeType;
    using typename BaseMesh::RealType;

    static constexpr Cubism::MeshClass Class = Cubism::MeshClass::Uniform;

private:
    using BaseMesh::crange_;
    using BaseMesh::frange_;
    using BaseMesh::nrange_;
    using BaseMesh::range_;

public:
    /**
     * @brief Standard mesh constructor
     * @param start Lower left point of physical domain
     * @param end Upper right point of physical domain
     * @param cells Number of cells in mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     */
    StructuredUniform(const PointType &start,
                      const PointType &end,
                      const MultiIndex &cells,
                      const MeshIntegrity type)
        : BaseMesh(start, end, cells, type),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    /**
     * @brief Standard mesh constructor
     * @param end Upper right point of physical domain
     * @param cells Number of cells in mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     *
     * Physical origin starts at 0
     */
    StructuredUniform(const PointType &end,
                      const MultiIndex &cells,
                      const MeshIntegrity type)
        : BaseMesh(end, cells, type),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    /**
     * @brief Standard mesh constructor
     * @param gorigin Global domain origin
     * @param range Domain range spanned by this mesh
     * @param crange Cell range spanned by this mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     *
     * Used for MPI subdomains
     */
    StructuredUniform(const PointType &gorigin,
                      const RangeType &range,
                      const IndexRangeType &crange,
                      const MeshIntegrity type)
        : BaseMesh(gorigin, range, crange, type),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    /**
     * @brief Low-level mesh constructor
     * @param gorigin Global domain origin
     * @param range Domain range spanned by this mesh
     * @param crange Cell range spanned by this mesh
     * @param nrange Node range spanned by this mesh
     * @param frange Face range spanned by this mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     *
     * Used for grid topology classes and sub-meshes
     */
    StructuredUniform(const PointType &gorigin,
                      const RangeType &range,
                      const IndexRangeType &crange,
                      const IndexRangeType &nrange,
                      const std::vector<IndexRangeType> &frange,
                      const MeshIntegrity type)
        : BaseMesh(gorigin, range, crange, nrange, frange, type),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    StructuredUniform() = delete;
    StructuredUniform(const StructuredUniform &c) = default;
    StructuredUniform(StructuredUniform &&c) = default;
    StructuredUniform &operator=(const StructuredUniform &c) = delete;
    StructuredUniform &operator=(StructuredUniform &&c) = delete;
    ~StructuredUniform() = default;

    /**
     * @brief Get sub-mesh instance
     * @param start Lower left point of physical domain for sub-mesh
     * @param end Upper right point of physical domain for sub-mesh
     * @return New sub-mesh instance
     *
     * @rst
     * If ``start`` or ``end`` is outside of the physical range spanned by this
     * mesh, then they will be adjusted to the start and/or end points of the
     * range corresponding to this mesh, respectively.  The extracted discrete
     * mesh is always guaranteed to include the ``start`` and ``end`` points.
     * Therefore, ``getRange().getBegin()`` of the extracted sub-mesh may be
     * smaller in any component than ``start`` and vice versa
     * ``getRange().getEnd()`` may be larger in any component than ``end``.
     *
     * .. note:: This method is not part of the ``BaseMesh`` interface.
     * @endrst
     */
    virtual std::unique_ptr<StructuredUniform>
    getSubMesh(const PointType &start, const PointType &end) const
    {
        const IndexRangeType sub_crange = getSubCellRange_(start, end);
        const PointType sub_start = getCoords_(
            sub_crange.getBegin() - crange_.getBegin(), EntityType::Node, 0);
        const PointType sub_end = getCoords_(
            (sub_crange.getEnd() - crange_.getBegin()), EntityType::Node, 0);
        const RangeType sub_range(sub_start, sub_end);
        return std::unique_ptr<StructuredUniform>(
            new StructuredUniform(this->getGlobalOrigin(),
                                  sub_range,
                                  sub_crange,
                                  MeshIntegrity::SubMesh));
    }

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

    IndexRangeType getSubCellRange_(PointType start, PointType end) const
    {
        const PointType end_init = end;
        const PointType start_r = range_.getBegin();
        const PointType end_r = range_.getEnd();
        if (start > end) {
            throw std::runtime_error("StructuredUniform: Can not create "
                                     "sub-cell range for start > end");
        } else if (!(start_r <= end && start <= end_r)) {
            throw std::runtime_error(
                "StructuredUniform: Range spanned by start and end points is "
                "not intersecting mesh range");
        }

        for (size_t i = 0; i < DIM; ++i) {
            // strip to this mesh range
            start[i] = (start[i] < start_r[i]) ? start_r[i] : start[i];
            end[i] = (end[i] > end_r[i]) ? end_r[i] : end[i];

            // boundaries: this ensures the correct resulting cell index if
            // start (and consequently end) happen to lie in the same cell
            // adjacent to the boundary
            start[i] = (Cubism::myAbs(start[i] - end_r[i]) < mesh_spacing_[i])
                           ? start[i] - 1.0 * mesh_spacing_[i]
                           : start[i];

            end[i] = (Cubism::myAbs(end[i] - start_r[i]) < mesh_spacing_[i])
                         ? end[i] + 1.0 * mesh_spacing_[i]
                         : end[i];

            // round up (internal region only): round up to closest cell index
            // (this ensures that the requested physical region is contained in
            // the extracted discrete index space)
            const int idx = (end[i] - start_r[i]) / mesh_spacing_[i];
            const RealType ep1 = start_r[i] + (idx + 1) * mesh_spacing_[i];
            const RealType adiff = Cubism::myAbs(ep1 - end[i]);
            const RealType inner = Cubism::myAbs(end_r[i] - end[i]);
            end[i] = (inner > mesh_spacing_[i] && 0 < adiff &&
                      adiff < mesh_spacing_[i])
                         ? end[i] + 1.0 * mesh_spacing_[i]
                         : end[i];

            // ensure at least one cell thick: the extracted region must be at
            // least one cell thick.  This criterion is ensured by always
            // rounding up.
            const RealType diff_rel = Cubism::myAbs(end[i] - start[i]);
            end[i] = (diff_rel < mesh_spacing_[i])
                         ? end[i] + 1.0 * mesh_spacing_[i]
                         : end[i];

        }
        MultiIndex ds((start - range_.getBegin()) / mesh_spacing_);
        MultiIndex de((end - range_.getBegin()) / mesh_spacing_);
        return IndexRangeType(crange_.getBegin() + ds, crange_.getBegin() + de);
    }
};

template <typename TReal, size_t DIM>
constexpr Cubism::MeshClass StructuredUniform<TReal, DIM>::Class;

NAMESPACE_END(Mesh)
NAMESPACE_END(Cubism)

#endif /* STRUCTUREDUNIFORM_H_L4AIO9HT */
