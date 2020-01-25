// File       : StructuredBase.h
// Created    : Thu Jan 02 2020 09:22:10 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic mesh interface
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef STRUCTUREDBASE_H_QETR3VEW
#define STRUCTUREDBASE_H_QETR3VEW

#include "Common.h"
#include "Core/Index.h"
#include "Core/Range.h"

#include <cassert>
#include <iterator>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Mesh)

/**
 * @brief Structured mesh base class
 * @tparam TReal Float type for mesh entities
 * @tparam DIM Mesh dimension
 *
 * Defines the (pure virtual) interface for a structured mesh.
 * */
template <typename TReal, size_t DIM>
class StructuredBase
{
public:
    /** @brief Range that spans the rectangular domain */
    using RangeType = Core::Range<TReal, DIM>;
    /** @brief Float type to describe mesh entities */
    using RealType = typename RangeType::DataType;
    /** @brief Point type in mesh */
    using PointType = typename RangeType::PointType;
    /** @brief Index range used to address discrete mesh entities */
    using IndexRangeType = Core::IndexRange<DIM>;
    /** @brief Multi-dimensional index */
    using MultiIndex = typename IndexRangeType::MultiIndex;
    /** @brief Mesh entity type */
    using EntityType = Cubism::EntityType;
    /** @brief Mesh integrity */
    using MeshIntegrity = Cubism::MeshIntegrity;

    static constexpr size_t Dim = DIM;

private:
    /** @brief Defines iterator classes for a ``Cubism::EntityType`` */
    class EntityIterator
    {
    public:
        EntityIterator(const EntityType t = EntityType::Cell,
                       const size_t d = 0,
                       const IndexRangeType &r = IndexRangeType())
            : entity_(t), dir_(d), range_(r)
        {
        }
        EntityIterator(const EntityIterator &c) = default;
        ~EntityIterator() = default;
        EntityIterator &operator=(const EntityIterator &c) = default;

        using iterator = Core::EntityIterator<DIM>;
        iterator begin() noexcept { return iterator(entity_, dir_, range_, 0); }
        iterator begin() const noexcept
        {
            return iterator(entity_, dir_, range_, 0);
        }
        iterator end() noexcept
        {
            return iterator(entity_, dir_, range_, range_.size());
        }
        iterator end() const noexcept
        {
            return iterator(entity_, dir_, range_, range_.size());
        }
        // TODO: [fabianw@mavt.ethz.ch; 2020-01-03] reverse iterator?

    private:
        const EntityType entity_;
        const size_t dir_;
        const IndexRangeType range_;
    };

    /**
     * @brief Iterator entity container for entity in ``Cubism::EntityType``
     * (for ``operator[]``)
     *
     * @rst
     * Faces are treated per dimension, which can be accessed with
     * ``Cubism::Dir::X`` for example.
     * @endrst
     */
    class Entity
    {
    public:
        Entity(const EntityType t, const IndexRangeType r = IndexRangeType())
            : ranges_(1, r), entities_(1, t)
        {
        }
        Entity(const EntityType t, const std::vector<IndexRangeType> &r)
            : ranges_(r), entities_(r.size(), t)
        {
        }
        Entity() = default;
        Entity(const Entity &c) = default;
        ~Entity() = default;
        Entity &operator=(const Entity &c) = default;

        using iterator = Core::EntityIterator<DIM>;
        iterator begin() noexcept
        {
            assert(ranges_.size() > 0);
            return iterator(entities_[0], 0, ranges_[0], 0);
        }
        iterator begin() const noexcept
        {
            assert(ranges_.size() > 0);
            return iterator(entities_[0], 0, ranges_[0], 0);
        }
        iterator end() noexcept
        {
            assert(ranges_.size() > 0);
            return iterator(entities_[0], 0, ranges_[0], ranges_[0].size());
        }
        iterator end() const noexcept
        {
            assert(ranges_.size() > 0);
            return iterator(entities_[0], 0, ranges_[0], ranges_[0].size());
        }
        // TODO: [fabianw@mavt.ethz.ch; 2020-01-03] reverse iterator?

        EntityIterator operator[](const size_t dir) const
        {
            assert(ranges_.size() > 0);
            assert(dir < DIM);
            return EntityIterator(entities_[dir], dir, ranges_[dir]);
        }

        template <typename Dir>
        EntityIterator operator[](const Dir d) const
        {
            return this->operator[](static_cast<size_t>(d));
        }

    private:
        std::vector<IndexRangeType> ranges_;
        std::vector<EntityType> entities_;
    };

public:
    /**
     * @brief Standard mesh constructor
     * @param start Lower left point of physical domain
     * @param end Upper right point of physical domain
     * @param cells Number of cells in mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     */
    StructuredBase(const PointType &start,
                   const PointType &end,
                   const MultiIndex &cells,
                   const MeshIntegrity type)
        : type_(type), range_(start, end), global_origin_(range_.getBegin()),
          crange_(cells), nrange_(crange_.getBegin(), crange_.getEnd() + 1)
    {
        initFaceRange_(crange_);
    }

    /**
     * @brief Standard mesh constructor
     * @param end Upper right point of physical domain
     * @param cells Number of cells in mesh
     * @param type Mesh integrity type (full mesh or sub-mesh)
     *
     * Physical origin starts at 0
     */
    StructuredBase(const PointType &end,
                   const MultiIndex &cells,
                   const MeshIntegrity type)
        : type_(type), range_(end), global_origin_(range_.getBegin()),
          crange_(cells), nrange_(crange_.getBegin(), crange_.getEnd() + 1)
    {
        initFaceRange_(crange_);
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
    StructuredBase(const PointType &gorigin,
                   const RangeType &range,
                   const IndexRangeType &crange,
                   const MeshIntegrity type)
        : type_(type), range_(range), global_origin_(gorigin), crange_(crange),
          nrange_(crange_.getBegin(), crange_.getEnd() + 1)
    {
        initFaceRange_(crange_);
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
    StructuredBase(const PointType &gorigin,
                   const RangeType &range,
                   const IndexRangeType &crange,
                   const IndexRangeType &nrange,
                   const std::vector<IndexRangeType> &frange,
                   const MeshIntegrity type)
        : type_(type), range_(range), global_origin_(gorigin), crange_(crange),
          nrange_(nrange), frange_(frange)
    {
    }

    StructuredBase() = delete;
    StructuredBase(StructuredBase &&c) = default;
    virtual ~StructuredBase() = default;

    StructuredBase(const StructuredBase &c)
        : type_(c.type_), range_(c.range_), global_origin_(c.global_origin_),
          crange_(c.crange_), nrange_(c.nrange_), frange_(c.frange_)
    {
    }

    StructuredBase &operator=(const StructuredBase &c) = delete;
    StructuredBase &operator=(StructuredBase &&c) = delete;

    /** @brief Mesh entity iterator */
    using iterator = typename EntityIterator::iterator;

    /**
     * @brief Get iterator for a mesh entity
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Entity iterator class
     */
    EntityIterator getIterator(const EntityType t, const size_t d = 0)
    {
        if (t == EntityType::Cell) {
            return EntityIterator(t, d, crange_);
        } else if (t == EntityType::Node) {
            return EntityIterator(t, d, nrange_);
        } else if (t == EntityType::Face) {
            return EntityIterator(t, d, frange_[d]);
        } else {
            throw std::runtime_error(
                "StructuredBase::getIterator: Unknown entity type t");
        }
        return EntityIterator();
    }

    /**
     * @brief Get iterator for a mesh entity
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Entity iterator class
     *
     * @rst
     * Can be used with ``Cubism::Dir::X`` for example.
     * @endrst
     */
    template <typename Dir>
    EntityIterator getIterator(const EntityType t, const Dir d)
    {
        return getIterator(t, static_cast<size_t>(d));
    }

    /**
     * @brief Iterable mesh entities
     * @param t Entity type
     * @return Entity iterator container
     *
     * @rst
     * Assume ``m`` is a mesh instance.  The access operator can be used in the
     * following manner:
     *
     * .. code-block:: cpp
     *
     *    for (auto c : m[EntityType::Cell]) {} // c is a cell index
     *
     *    auto faces = m[EntityType::Face]; // iterators for faces
     *    for (auto fx : faces) {}         // fx is a X-face index
     *    for (auto fx : faces[Dir::X]) {} // fx is a X-face index
     *    for (auto fy : faces[Dir::Y]) {} // fy is a Y-face index
     *    for (auto fz : faces[Dir::Z]) {} // fz is a Z-face index
     * @endrst
     */
    Entity operator[](const EntityType t) const
    {
        if (t == EntityType::Cell) {
            return Entity(t, crange_);
        } else if (t == EntityType::Node) {
            return Entity(t, nrange_);
        } else if (t == EntityType::Face) {
            return Entity(t, frange_);
        } else {
            throw std::runtime_error(
                "StructuredBase::operator[]: Unknown entity type t");
        }
        return Entity();
    }

    /**
     * @brief Total size of mesh entities
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Number of entities
     */
    size_t size(const EntityType t, const size_t d = 0) const
    {
        assert(d < frange_.size());
        if (t == EntityType::Cell) {
            return crange_.size();
        } else if (t == EntityType::Node) {
            return nrange_.size();
        } else if (t == EntityType::Face) {
            return frange_[d].size();
        } else {
            throw std::runtime_error(
                "StructuredBase::size: Unknown entity type t");
        }
        return 0;
    }

    /**
     * @brief Total size of mesh entities
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Number of entities
     */
    template <typename Dir>
    size_t size(const EntityType t, const Dir d) const
    {
        return size(t, static_cast<size_t>(d));
    }

    /**
     * @brief Size of mesh entities
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Number of entities along all directions
     */
    MultiIndex getSize(const EntityType t, const size_t d = 0) const
    {
        assert(d < frange_.size());
        if (t == EntityType::Cell) {
            return crange_.getExtent();
        } else if (t == EntityType::Node) {
            return nrange_.getExtent();
        } else if (t == EntityType::Face) {
            return frange_[d].getExtent();
        } else {
            throw std::runtime_error(
                "StructuredBase::getSize: Unknown entity type t");
        }
        return MultiIndex();
    }

    /**
     * @brief Size of mesh entities
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Number of entities along all directions
     */
    template <typename Dir>
    MultiIndex getSize(const EntityType t, const Dir d) const
    {
        return getSize(t, static_cast<size_t>(d));
    }

    /**
     * @brief Index range for an entity type
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Index range
     */
    IndexRangeType getIndexRange(const EntityType t, const size_t d = 0) const
    {
        assert(d < frange_.size());
        if (t == EntityType::Cell) {
            return crange_;
        } else if (t == EntityType::Node) {
            return nrange_;
        } else if (t == EntityType::Face) {
            return frange_[d];
        } else {
            throw std::runtime_error(
                "StructuredBase::getIndexRange: Unknown entity type t");
        }
        return IndexRangeType();
    }

    /**
     * @brief Index range for an entity type
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Index range
     */
    template <typename Dir>
    IndexRangeType getIndexRange(const EntityType t, const Dir d) const
    {
        return getIndexRange(t, static_cast<size_t>(d));
    }

    /**
     * @brief Mesh extent
     * @return Extent of mesh in all directions
     */
    PointType getExtent() const { return range_.getExtent(); }
    /**
     * @brief Mesh volume
     * @return Total mesh volume
     */
    RealType getVolume() const { return range_.getVolume(); }
    /**
     * @brief Get mesh origin
     * @return Local origin of the mesh
     */
    PointType getOrigin() const { return range_.getBegin(); }
    /**
     * @brief Get mesh origin
     * @return Global origin of the mesh
     */
    PointType getGlobalOrigin() const { return global_origin_; }
    /**
     * @brief Get mesh range
     * @return Mesh range
     */
    RangeType getRange() const { return range_; }

    /**
     * @brief Test if this is a sub-mesh
     * @return True if this is a sub-mesh of some full mesh
     */
    bool isSubMesh() const { return type_ == MeshIntegrity::SubMesh; }

    /**
     * @brief Convert local flat index to multi-dimensional index
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Local multi-dimensional index
     */
    MultiIndex
    getMultiIndex(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getMultiIndex_(i, t, d);
    }

    /**
     * @brief Convert local flat index to multi-dimensional index
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Local multi-dimensional index
     */
    template <typename Dir>
    MultiIndex
    getMultiIndex(const size_t i, const EntityType t, const Dir d) const
    {
        return getMultiIndex_(i, t, static_cast<size_t>(d));
    }

    /**
     * @brief Convert local flat index to global multi-dimensional index
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global multi-dimensional index
     */
    MultiIndex
    getGlobalIndex(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getGlobalIndex(getMultiIndex_(i, t, d), t, d);
    }

    /**
     * @brief Convert local flat index to global multi-dimensional index
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global multi-dimensional index
     */
    template <typename Dir>
    MultiIndex
    getGlobalIndex(const size_t i, const EntityType t, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getGlobalIndex(getMultiIndex_(i, t, di), t, di);
    }

    /**
     * @brief Convert local multi-dimensional index to global multi-dimensional
     * index
     * @param p Local multi-dimensional index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global multi-dimensional index
     */
    MultiIndex getGlobalIndex(const MultiIndex &p,
                              const EntityType t,
                              const size_t d = 0) const
    {
        assert(d < frange_.size());
        if (t == EntityType::Cell) {
            return crange_.getBegin() + p;
        } else if (t == EntityType::Node) {
            return nrange_.getBegin() + p;
        } else if (t == EntityType::Face) {
            return frange_[d].getBegin() + p;
        } else {
            throw std::runtime_error(
                "StructuredBase::getGlobalIndex: Unknown entity type t");
        }
        return p;
    }

    /**
     * @brief Convert local multi-dimensional index to global multi-dimensional
     * index
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global multi-dimensional index
     */
    template <typename Dir>
    MultiIndex
    getGlobalIndex(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return getGlobalIndex(p, t, static_cast<size_t>(d));
    }

    /**
     * @brief Convert local flat index to global coordinates
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global coordinates
     */
    PointType getGlobalCoords(const size_t i,
                              const EntityType t,
                              const size_t d = 0) const
    {
        assert(d < frange_.size());
        return global_origin_ + getCoords(i, t, d);
    }

    /**
     * @brief Convert local flat index to global coordinates
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global coordinates
     */
    template <typename Dir>
    PointType
    getGlobalCoords(const size_t i, const EntityType t, const Dir d) const
    {
        return global_origin_ + getCoords(i, t, d);
    }

    /**
     * @brief Convert local multi-dimensional index to global coordinates
     * @param p Local multi-dimensional index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global coordinates
     */
    PointType getGlobalCoords(const MultiIndex &p,
                              const EntityType t,
                              const size_t d = 0) const
    {
        assert(d < frange_.size());
        return global_origin_ + getCoords(p, t, d);
    }

    /**
     * @brief Convert local multi-dimensional index to global coordinates
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Global coordinates
     */
    template <typename Dir>
    PointType
    getGlobalCoords(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return global_origin_ + getCoords(p, t, d);
    }

    /**
     * @brief Convert entity iterator to global coordinates
     * @param it Entity iterator
     * @return Global coordinates
     */
    PointType getGlobalCoords(const iterator &it) const
    {
        return global_origin_ + getCoords(it);
    }

    /**
     * @brief Convert local flat index to local coordinates
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Local coordinates
     */
    PointType
    getCoords(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getCoords_(getMultiIndex_(i, t, d), t, d);
    }

    /**
     * @brief Convert local flat index to local coordinates
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Local coordinates
     */
    template <typename Dir>
    PointType getCoords(const size_t i, const EntityType t, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getCoords_(getMultiIndex_(i, t, di), t, di);
    }

    /**
     * @brief Convert local multi-dimensional index to local coordinates
     * @param p Local multi-dimensional index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Local coordinates
     */
    PointType
    getCoords(const MultiIndex &p, const EntityType t, const size_t d = 0) const
    {
        assert(d < frange_.size());
        return getCoords_(p, t, d);
    }

    /**
     * @brief Convert local multi-dimensional index to local coordinates
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional index
     * @param t Entity type
     * @param d Direction indicator (for Face entities only)
     * @return Local coordinates
     */
    template <typename Dir>
    PointType
    getCoords(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return getCoords_(p, t, static_cast<size_t>(d));
    }

    /**
     * @brief Convert entity iterator to local coordinates
     * @param it Entity iterator
     * @return Local coordinates
     */
    PointType getCoords(const iterator &it) const
    {
        return getCoords_(*it, it.getEntity(), it.getDirection());
    }

    /**
     * @brief Get cell volume
     * @param i Local flat cell index
     * @return Cell volume
     */
    RealType getCellVolume(const size_t i) const
    {
        return getCellVolume_(getMultiIndex_(i, EntityType::Cell));
    }

    /**
     * @brief Get cell volume
     * @param p Local multi-dimensional cell index
     * @return Cell volume
     */
    RealType getCellVolume(const MultiIndex &p) const
    {
        return getCellVolume_(p);
    }

    /**
     * @brief Get cell volume
     * @param it Cell entity iterator
     * @return Cell volume
     */
    RealType getCellVolume(const iterator &it) const
    {
        return getCellVolume_(*it);
    }

    /**
     * @brief Get cell size
     * @param i Local flat cell index
     * @return Cell size
     */
    PointType getCellSize(const size_t i) const
    {
        return getCellSize_(getMultiIndex_(i, EntityType::Cell));
    }

    /**
     * @brief Get cell size
     * @param p Local multi-dimensional cell index
     * @return Cell size
     */
    PointType getCellSize(const MultiIndex &p) const { return getCellSize_(p); }

    /**
     * @brief Get cell size
     * @param it Cell entity iterator
     * @return Cell size
     */
    PointType getCellSize(const iterator &it) const
    {
        return getCellSize_(*it);
    }

    /**
     * @brief Get surface vector for a face
     * @param i Local flat face index
     * @param ci Local flat index of reference cell
     * @param d Direction of face
     * @return Face surface vector
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    PointType getSurface(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d);
    }

    /**
     * @brief Get surface vector for a face
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat face index
     * @param ci Local flat index of reference cell
     * @param d Direction of face
     * @return Face surface vector
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    template <typename Dir>
    PointType getSurface(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di);
    }

    /**
     * @brief Get surface vector for a face
     * @param p Local multi-dimensional face index
     * @param ci Local multi-dimensional index of reference cell
     * @param d Direction of face
     * @return Face surface vector
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    PointType
    getSurface(const MultiIndex &p, const MultiIndex &ci, const size_t d) const
    {
        assert(d < frange_.size());
        return getSurface_(p, ci, d);
    }

    /**
     * @brief Get surface vector for a face
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional face index
     * @param ci Local multi-dimensional index of reference cell
     * @param d Direction of face
     * @return Face surface vector
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    template <typename Dir>
    PointType
    getSurface(const MultiIndex &p, const MultiIndex &ci, const Dir d) const
    {
        return getSurface_(p, ci, static_cast<size_t>(d));
    }

    // FIXME: [fabianw@mavt.ethz.ch; 2020-01-04] I don't think this interface is
    // very useful.  Removal candidate
    PointType getSurface(const iterator &it, const MultiIndex &ci) const
    {
        return getSurface_(*it, ci, it.getDirection());
    }

    /**
     * @brief Get surface area for a face
     * @param i Local flat face index
     * @param ci Local flat index of reference cell
     * @param d Direction of face
     * @return Face surface area
     */
    RealType
    getSurfaceArea(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d)
            .norm();
    }

    /**
     * @brief Get surface area for a face
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat face index
     * @param ci Local flat index of reference cell
     * @param d Direction of face
     * @return Face surface area
     */
    template <typename Dir>
    RealType getSurfaceArea(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di)
            .norm();
    }

    /**
     * @brief Get surface area for a face
     * @param p Local multi-dimensional face index
     * @param ci Local multi-dimensional index of reference cell
     * @param d Direction of face
     * @return Face surface area
     */
    RealType getSurfaceArea(const MultiIndex &p,
                            const MultiIndex &ci,
                            const size_t d) const
    {
        assert(d < frange_.size());
        return getSurface_(p, ci, d).norm();
    }

    /**
     * @brief Get surface area for a face
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional face index
     * @param ci Local multi-dimensional index of reference cell
     * @param d Direction of face
     * @return Face surface area
     */
    template <typename Dir>
    RealType
    getSurfaceArea(const MultiIndex &p, const MultiIndex &ci, const Dir d) const
    {
        return getSurface_(p, ci, static_cast<size_t>(d)).norm();
    }

    // FIXME: [fabianw@mavt.ethz.ch; 2020-01-04] I don't think this interface is
    // very useful.  Removal candidate
    RealType getSurfaceArea(const iterator &it, const MultiIndex &ci) const
    {
        return getSurface_(*it, ci, it.getDirection()).norm();
    }

    /**
     * @brief Get surface normal for a face
     * @param i Local flat face index
     * @param ci Local flat index of reference cell
     * @param d Direction of face
     * @return Face surface normal
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    PointType
    getSurfaceNormal(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d)
            .unit();
    }

    /**
     * @brief Get surface normal for a face
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param i Local flat face index
     * @param ci Local flat index of reference cell
     * @param d Direction of face
     * @return Face surface normal
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    template <typename Dir>
    PointType
    getSurfaceNormal(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di)
            .unit();
    }

    /**
     * @brief Get surface normal for a face
     * @param p Local multi-dimensional face index
     * @param ci Local multi-dimensional index of reference cell
     * @param d Direction of face
     * @return Face surface normal
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    PointType getSurfaceNormal(const MultiIndex &p,
                               const MultiIndex &ci,
                               const size_t d) const
    {
        assert(d < frange_.size());
        return getSurface_(p, ci, d).unit();
    }

    /**
     * @brief Get surface normal for a face
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param p Local multi-dimensional face index
     * @param ci Local multi-dimensional index of reference cell
     * @param d Direction of face
     * @return Face surface normal
     *
     * @rst
     * The vector points *outward* of the specified reference cell ``ci``.
     * @endrst
     */
    template <typename Dir>
    PointType getSurfaceNormal(const MultiIndex &p,
                               const MultiIndex ci,
                               const Dir d) const
    {
        return getSurface_(p, ci, static_cast<size_t>(d)).unit();
    }

    // FIXME: [fabianw@mavt.ethz.ch; 2020-01-04] I don't think this interface is
    // very useful.  Removal candidate
    PointType getSurfaceNormal(const iterator &it, const MultiIndex &ci) const
    {
        return getSurface_(*it, ci, it.getDirection()).unit();
    }

protected:
    const MeshIntegrity type_;      // mesh integrity type (full or sub-mesh)
    const RangeType range_;         // range of mesh domain in physical space
    const PointType global_origin_; // origin in global mesh
    const IndexRangeType crange_;   // index range spanned by cells
    const IndexRangeType nrange_;   // index range spanned by nodes
    std::vector<IndexRangeType> frange_; // index ranges spanned by faces

    /**
     * @brief Pure virtual method to obtain local coordinates
     * @param p Local multi-dimensional index
     * @param t Cubism::EntityType
     * @param dir Entity direction (for faces only)
     */
    virtual PointType getCoords_(const MultiIndex &p,
                                 const EntityType t,
                                 const size_t dir) const = 0;

    /**
     * @brief Pure virtual method to obtain cell volume
     * @param p Local multi-dimensional index
     */
    virtual RealType getCellVolume_(const MultiIndex &p) const = 0;

    /**
     * @brief Pure virtual method to obtain cell size
     * @param p Local multi-dimensional index
     */
    virtual PointType getCellSize_(const MultiIndex &p) const = 0;

    /**
     * @brief Pure virtual method to obtain surface vector
     * @param p Local multi-dimensional index
     * @param ci Local multi-dimensional index of reference cell (determines
     *           orientation)
     * @param dir Entity direction (for faces only)
     */
    virtual PointType getSurface_(const MultiIndex &p,
                                  const MultiIndex &ci,
                                  const size_t dir) const = 0;

private:
    void initFaceRange_(const IndexRangeType &r)
    {
        for (size_t i = 0; i < DIM; ++i) {
            frange_.push_back(IndexRangeType(
                r.getBegin(), r.getEnd() + MultiIndex::getUnitVector(i)));
        }
    }

    MultiIndex
    getMultiIndex_(const size_t i, const EntityType t, const size_t d = 0) const
    {
        assert(d < frange_.size());
        if (t == EntityType::Cell) {
            return crange_.getMultiIndex(i);
        } else if (t == EntityType::Node) {
            return nrange_.getMultiIndex(i);
        } else if (t == EntityType::Face) {
            return frange_[d].getMultiIndex(i);
        } else {
            throw std::runtime_error(
                "StructuredBase::getMultiIndex_: Unknown entity type t");
        }
        return MultiIndex();
    }
};

template <typename TReal, size_t DIM>
constexpr size_t StructuredBase<TReal, DIM>::Dim;

NAMESPACE_END(Mesh)
NAMESPACE_END(Cubism)

#endif /* STRUCTUREDBASE_H_QETR3VEW */
