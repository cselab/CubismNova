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

/**
 * @addtogroup Mesh
 * @{
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

template <typename TReal, size_t DIM>
constexpr Cubism::MeshClass StructuredUniform<TReal, DIM>::Class;

/**  @} */
NAMESPACE_END(Mesh)
NAMESPACE_END(Cubism)

#endif /* STRUCTUREDUNIFORM_H_L4AIO9HT */
