// File       : StructuredUniform.h
// Created    : Sat Jan 04 2020 11:17:40 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Structured uniform mesh
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef STRUCTUREDUNIFORM_H_L4AIO9HT
#define STRUCTUREDUNIFORM_H_L4AIO9HT

#include "Cubism/Math.h"
#include "Cubism/Mesh/StructuredBase.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Mesh)

/**
 * @ingroup Mesh
 * @brief Structured uniform mesh
 * @tparam TReal Float type for mesh entities
 * @tparam DIM Mesh dimension
 * */
template <typename TReal, size_t DIM = CUBISM_DIMENSION>
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
    using BaseMesh::global_range_;
    using BaseMesh::nrange_;
    using BaseMesh::range_;

public:
    /**
     * @brief Standard mesh constructor
     * @param begin Lower left point of physical domain
     * @param end Upper right point of physical domain
     * @param cells Number of cells in mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     */
    StructuredUniform(const PointType &begin,
                      const PointType &end,
                      const MultiIndex &cells,
                      const MeshIntegrity type)
        : BaseMesh(begin, end, cells, type),
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
     * Physical begin at 0
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
     * @param cells Number of cells in mesh
     *
     * Physical domain [0, 1], always a full mesh type
     */
    StructuredUniform(const MultiIndex &cells)
        : BaseMesh(cells),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    /**
     * @brief Standard mesh constructor
     * @param grange Global domain range spanned by this mesh
     * @param range Domain range spanned by this mesh
     * @param crange Cell range spanned by this mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     *
     * Used for MPI subdomains
     */
    StructuredUniform(const RangeType &grange,
                      const RangeType &range,
                      const IndexRangeType &crange,
                      const MeshIntegrity type)
        : BaseMesh(grange, range, crange, type),
          mesh_spacing_(range_.getExtent() / PointType(crange_.getExtent())),
          cell_volume_(mesh_spacing_.prod())
    {
    }

    /**
     * @brief Low-level mesh constructor
     * @param grange Global domain range spanned by this mesh
     * @param range Domain range spanned by this mesh
     * @param crange Cell range spanned by this mesh
     * @param nrange Node range spanned by this mesh
     * @param frange Face range spanned by this mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     *
     * Used for grid topology classes and sub-meshes
     */
    StructuredUniform(const RangeType &grange,
                      const RangeType &range,
                      const IndexRangeType &crange,
                      const IndexRangeType &nrange,
                      const std::vector<IndexRangeType> &frange,
                      const MeshIntegrity type)
        : BaseMesh(grange, range, crange, nrange, frange, type),
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
     * @param range Index range of new mesh
     * @param entity Entity type corresponding to ``range``
     * @param d Direction indicator (for Face entities only)
     * @return New sub-mesh instance
     *
     * @rst
     * An empty mesh is returned if there is no common intersection.  If
     * ``range`` is larger it will be clipped to the boundaries of this mesh.
     *
     * .. note:: This method is not part of the ``BaseMesh`` interface.
     * @endrst
     */
    virtual std::unique_ptr<StructuredUniform>
    getSubMesh(const IndexRangeType &range,
               const EntityType entity,
               const size_t d = 0) const
    {
        IndexRangeType common =
            this->getIndexRange(entity, d).getIntersection(range);
        const auto null = common.getNullSpace();
        MultiIndex cend = common.getEnd();
        if (null.size() != DIM) {
            if (entity == Cubism::EntityType::Node) {
                cend -= 1;
                for (const auto &i : null) {
                    cend[i] += 1; // maintain degenerate space
                }
            } else if (entity == Cubism::EntityType::Face) {
                if (std::find(null.begin(), null.end(), d) == null.end()) {
                    cend -= MultiIndex::getUnitVector(d);
                }
            }
        }
        common.setEnd(cend); // cell range now
        const PointType sub_begin = getCoords_(
            common.getBegin() - crange_.getBegin(), EntityType::Node, 0);
        if (cend != crange_.getEnd()) {
            cend += 1; // top right node index
            for (const auto &i : null) {
                cend[i] -= 1; // maintain degenerate space
            }
        }
        const PointType sub_end =
            getCoords_(cend - crange_.getBegin(), EntityType::Node, 0);
        const RangeType sub_range(sub_begin, sub_end);
        return std::unique_ptr<StructuredUniform>(new StructuredUniform(
            this->getGlobalRange(), sub_range, common, MeshIntegrity::SubMesh));
    }

    /**
     * @brief Get sub-mesh instance
     * @param begin Lower left point of physical domain for sub-mesh
     * @param end Upper right point of physical domain for sub-mesh
     * @return New sub-mesh instance
     *
     * @rst
     * If ``begin`` or ``end`` is outside of the physical range spanned by this
     * mesh, then they will be adjusted to the begin and/or end points of the
     * range corresponding to this mesh, respectively.  The extracted discrete
     * mesh is always guaranteed to include the ``begin`` and ``end`` points.
     * Therefore, ``getRange().getBegin()`` of the extracted sub-mesh may be
     * smaller in any component than ``begin`` and vice versa
     * ``getRange().getEnd()`` may be larger in any component than ``end``.
     *
     * .. note:: This method is not part of the ``BaseMesh`` interface.
     * @endrst
     */
    virtual std::unique_ptr<StructuredUniform>
    getSubMesh(const PointType &begin, const PointType &end) const
    {
        const IndexRangeType sub_crange = getSubCellRange_(begin, end);
        const PointType sub_begin = getCoords_(
            sub_crange.getBegin() - crange_.getBegin(), EntityType::Node, 0);
        const PointType sub_end = getCoords_(
            sub_crange.getEnd() - crange_.getBegin(), EntityType::Node, 0);
        const RangeType sub_range(sub_begin, sub_end);
        return std::unique_ptr<StructuredUniform>(
            new StructuredUniform(this->getGlobalRange(),
                                  sub_range,
                                  sub_crange,
                                  MeshIntegrity::SubMesh));
    }

    // non-virtual overloads for coordinate lookup:
    /**
     * @brief Get global cell coordinates
     * @param p Local multi-dimensional cell index
     * @return Global cell coordinates
     *
     * This is a non-virtual method. Prefer this method for excessive coordinate
     * lookup in loops.
     */
    PointType getGlobalCoordsCell(const MultiIndex &p) const
    {
        return global_range_.getBegin() + getCoordsCell_(p);
    }
    /**
     * @brief Get global node coordinates
     * @param p Local multi-dimensional node index
     * @return Global node coordinates
     *
     * This is a non-virtual method. Prefer this method for excessive coordinate
     * lookup in loops.
     */
    PointType getGlobalCoordsNode(const MultiIndex &p) const
    {
        return global_range_.getBegin() + getCoordsNode_(p);
    }
    /**
     * @brief Get global face coordinates
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional face index
     * @param dir Face direction identifier
     * @return Global face coordinates
     *
     * This is a non-virtual method. Prefer this method for excessive coordinate
     * lookup in loops.
     */
    template <typename Dir = size_t>
    PointType getGlobalCoordsFace(const MultiIndex &p, const Dir dir) const
    {
        return global_range_.getBegin() + getCoordsFace_(p, dir);
    }
    /**
     * @brief Get local cell coordinates
     * @param p Local multi-dimensional cell index
     * @return Local cell coordinates
     *
     * This is a non-virtual method. Prefer this method for excessive coordinate
     * lookup in loops.
     */
    PointType getCoordsCell(const MultiIndex &p) const
    {
        return getCoordsCell_(p);
    }
    /**
     * @brief Get local node coordinates
     * @param p Local multi-dimensional node index
     * @return Local node coordinates
     *
     * This is a non-virtual method. Prefer this method for excessive coordinate
     * lookup in loops.
     */
    PointType getCoordsNode(const MultiIndex &p) const
    {
        return getCoordsNode_(p);
    }
    /**
     * @brief Get local face coordinates
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional face index
     * @param dir Face direction identifier
     * @return Local face coordinates
     *
     * This is a non-virtual method. Prefer this method for excessive coordinate
     * lookup in loops.
     */
    template <typename Dir = size_t>
    PointType getCoordsFace(const MultiIndex &p, const Dir dir) const
    {
        return getCoordsFace_(p, dir);
    }

protected:
    PointType getCoords_(const MultiIndex &p,
                         const EntityType t,
                         const size_t dir) const override
    {
        if (t == EntityType::Cell) {
            return getCoordsCell_(p);
        } else if (t == EntityType::Node) {
            return getCoordsNode_(p);
        } else if (t == EntityType::Face) {
            return getCoordsFace_(p, dir);
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

    PointType getCoordsCell_(const MultiIndex &p) const
    {
        PointType c = PointType(p) + 0.5;
        return range_.getBegin() + c * mesh_spacing_;
    }

    PointType getCoordsNode_(const MultiIndex &p) const
    {
        PointType c(p);
        return range_.getBegin() + c * mesh_spacing_;
    }

    template <typename Dir = size_t>
    PointType getCoordsFace_(const MultiIndex &p, const Dir dir) const
    {
        PointType c =
            PointType(p) +
            0.5 * (1.0 - PointType::getUnitVector(static_cast<size_t>(dir)));
        return range_.getBegin() + c * mesh_spacing_;
    }

    IndexRangeType getSubCellRange_(PointType begin, PointType end) const
    {
        const PointType end_init = end;
        const PointType begin_r = range_.getBegin();
        const PointType end_r = range_.getEnd();
        if (begin > end) {
            throw std::runtime_error("StructuredUniform: Can not create "
                                     "sub-cell range for begin > end");
        } else if (!(begin_r <= end && begin <= end_r)) {
            return IndexRangeType();
        }

        for (size_t i = 0; i < DIM; ++i) {
            // strip to this mesh range
            begin[i] = (begin[i] < begin_r[i]) ? begin_r[i] : begin[i];
            end[i] = (end[i] > end_r[i]) ? end_r[i] : end[i];

            // boundaries: this ensures the correct cell index if begin (and
            // consequently end) happen to lie in the same cell adjacent to the
            // boundary or on the boundary up to machine precision.
            bool boundary = false;
            const RealType badiff = Cubism::myAbs(begin[i] - end_r[i]);
            if (badiff < 2.0 * std::numeric_limits<RealType>::epsilon()) {
                // avoid floating point issues
                begin[i] -= 0.5 * mesh_spacing_[i];
                end[i] += 0.5 * mesh_spacing_[i];
                boundary = true;
            } else if (0 < badiff && badiff < mesh_spacing_[i]) {
                // somewhere in the cell adjacent to the boundary
                end[i] += 1.0 * mesh_spacing_[i];
                boundary = true;
            }

            // round up (internal region only): round up to closest cell index
            // (this ensures that the requested physical region is contained in
            // the extracted discrete index space)
            const int idx = (end[i] - begin_r[i]) / mesh_spacing_[i];
            const RealType ep1 = begin_r[i] + (idx + 1) * mesh_spacing_[i];
            const RealType adiff = Cubism::myAbs(ep1 - end[i]);
            const RealType inner = Cubism::myAbs(end_r[i] - end[i]);
            if (!boundary && inner > mesh_spacing_[i] && 0 < adiff &&
                adiff < mesh_spacing_[i]) {
                end[i] += 1.0 * mesh_spacing_[i];
            }

            // ensure at least one cell thick: the extracted region must be at
            // least one cell thick.  This criterion is ensured by always
            // rounding up.
            const RealType diff_rel0 = Cubism::myAbs(end[i] - begin[i]);
            const RealType diff_rel1 =
                Cubism::myAbs(diff_rel0 - mesh_spacing_[i]);
            if (!boundary &&
                (diff_rel1 < 2.0 * std::numeric_limits<RealType>::epsilon())) {
                end[i] += 2.0 * std::numeric_limits<RealType>::epsilon();
            } else if (!boundary && (diff_rel0 < mesh_spacing_[i])) {
                end[i] += 1.0 * mesh_spacing_[i];
            }
        }
        MultiIndex ds((begin - range_.getBegin()) / mesh_spacing_);
        MultiIndex de((end - range_.getBegin()) / mesh_spacing_);
        return IndexRangeType(crange_.getBegin() + ds, crange_.getBegin() + de);
    }
};

template <typename TReal, size_t DIM>
constexpr Cubism::MeshClass StructuredUniform<TReal, DIM>::Class;

NAMESPACE_END(Mesh)
NAMESPACE_END(Cubism)

#endif /* STRUCTUREDUNIFORM_H_L4AIO9HT */
