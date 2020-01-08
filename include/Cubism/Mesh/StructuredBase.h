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

template <typename TReal, size_t DIM>
class StructuredBase
{
public:
    using RangeType = Core::Range<TReal, DIM>;
    using RealType = typename RangeType::DataType;
    using PointType = typename RangeType::PointType;
    using IndexRangeType = Core::IndexRange<DIM>;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using EntityType = Cubism::EntityType;

    enum class MeshHull { FullMesh = 0, SubMesh };
    enum class MeshClass { Uniform = 0, Stretched };

    static constexpr size_t Dim = DIM;

private:
    /// @brief Defines iterator classes for a Cubism::EntityType which are used
    ///        in range-based operators like for-loops
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

        using iterator = Core::MultiIndexIterator<DIM>;
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

    /// @brief Iterator entity for entity in Cubism::EntityType (for operator[])
    ///
    /// Faces are treated per dimension, which can be accessed with
    /// Cubism::Dir::X for example.
    class Entity
    {
    public:
        Entity(const EntityType t, const IndexRangeType r = IndexRangeType())
            : ranges_(1, r), entities_(1, t)
        {
        }
        Entity(const EntityType t, const std::vector<const IndexRangeType *> &r)
            : ranges_(), entities_(r.size(), t)
        {
            for (const auto p : r) {
                ranges_.push_back(*p);
            }
        }
        Entity() = default;
        Entity(const Entity &c) = default;
        ~Entity() = default;
        Entity &operator=(const Entity &c) = default;

        using iterator = Core::MultiIndexIterator<DIM>;
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
    /// @brief Standard mesh constructor
    ///
    /// @param start Lower left point of physical domain
    /// @param end Upper right point of physical domain
    /// @param cells Number of cells in mesh
    /// @param type Mesh hull type (full mesh or sub-mesh)
    /// @param cl Class of mesh (uniform, stretched)
    StructuredBase(const PointType &start,
                   const PointType &end,
                   const MultiIndex &cells,
                   const MeshHull type,
                   const MeshClass cl)
        : type_(type), class_(cl), range_(start, end),
          global_origin_(range_.getBegin()), crange_(cells),
          nrange_(crange_.getBegin(), crange_.getEnd() + 1),
          frange_(DIM, nullptr)
    {
        initFaceRange_(crange_);
    }

    /// @brief Standard mesh constructor (with physical origin at 0)
    ///
    /// @param end Upper right point of physical domain
    /// @param cells Number of cells in mesh
    /// @param type Mesh hull type (full mesh or sub-mesh)
    /// @param cl Class of mesh (uniform, stretched)
    StructuredBase(const PointType &end,
                   const MultiIndex &cells,
                   const MeshHull type,
                   const MeshClass cl)
        : type_(type), class_(cl), range_(end),
          global_origin_(range_.getBegin()), crange_(cells),
          nrange_(crange_.getBegin(), crange_.getEnd() + 1),
          frange_(DIM, nullptr)
    {
        initFaceRange_(crange_);
    }

    /// @brief Standard mesh constructor (useful for MPI subdomains)
    ///
    /// @param gorigin Global domain origin
    /// @param range Domain range spanned by this mesh
    /// @param crange Cell range spanned by this mesh
    /// @param type Mesh hull type (full mesh or sub-mesh)
    /// @param cl Class of mesh (uniform, stretched)
    StructuredBase(const PointType &gorigin,
                   const RangeType &range,
                   const IndexRangeType &crange,
                   const MeshHull type,
                   const MeshClass cl)
        : type_(type), class_(cl), range_(range), global_origin_(gorigin),
          crange_(crange), nrange_(crange_.getBegin(), crange_.getEnd() + 1),
          frange_(DIM, nullptr)
    {
        initFaceRange_(crange_);
    }

    /// @brief Low-level mesh constructor (useful for grid topology classes and
    ///        sub-meshes)
    ///
    /// @param gorigin Global domain origin
    /// @param range Domain range spanned by this mesh
    /// @param crange Cell range spanned by this mesh
    /// @param nrange Node range spanned by this mesh
    /// @param frange Face range spanned by this mesh
    /// @param type Mesh hull type (full mesh or sub-mesh)
    /// @param cl Class of mesh (uniform, stretched)
    StructuredBase(const PointType &gorigin,
                   const RangeType &range,
                   const IndexRangeType &crange,
                   const IndexRangeType &nrange,
                   const std::vector<IndexRangeType> &frange,
                   const MeshHull type,
                   const MeshClass cl)
        : type_(type), class_(cl), range_(range), global_origin_(gorigin),
          crange_(crange), nrange_(nrange), frange_(DIM, nullptr)
    {
        for (size_t i = 0; i < frange_.size(); ++i) {
            frange_[i] =
                new IndexRangeType(frange[i].getBegin(), frange[i].getEnd());
        }
    }

    StructuredBase() = delete;
    StructuredBase(StructuredBase &&c) = default;
    virtual ~StructuredBase() { disposeFaceRange_(); }

    StructuredBase(const StructuredBase &c)
        : type_(c.type_), class_(c.class_), range_(c.range_),
          global_origin_(c.global_origin_), crange_(c.crange_),
          nrange_(c.nrange_), frange_(DIM, nullptr)
    {
        copyFaceRange_(c);
    }

    StructuredBase &operator=(const StructuredBase &c) = delete;
    StructuredBase &operator=(StructuredBase &&c) = delete;

    // iterators
    using iterator = typename EntityIterator::iterator;

    /// @brief Returns iterator over Cubism::EntityType
    ///
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    EntityIterator getIterator(const EntityType t, const size_t d = 0)
    {
        if (t == EntityType::Cell) {
            return EntityIterator(t, d, crange_);
        } else if (t == EntityType::Node) {
            return EntityIterator(t, d, nrange_);
        } else if (t == EntityType::Face) {
            return EntityIterator(t, d, *frange_[d]);
        } else {
            throw std::runtime_error(
                "StructuredBase::getIterator: Unknown entity type t");
        }
        return EntityIterator();
    }

    /// @brief Returns iterator over Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    EntityIterator getIterator(const EntityType t, const Dir d)
    {
        return getIterator(t, static_cast<size_t>(d));
    }

    /// @brief Returns an entity iterator for Cubism::EntityType
    ///
    /// @param t Cubism::EntityType
    ///
    /// Examples:
    ///  - for (auto c : m[EntityType::Cell])
    ///  - for (auto f : m[EntityType::Face][Dir::X])
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

    /// @brief Returns the number of elements in the mesh for a
    ///        Cubism::EntityType
    ///
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    size_t size(const EntityType t, const size_t d = 0) const
    {
        if (t == EntityType::Cell) {
            return crange_.size();
        } else if (t == EntityType::Node) {
            return nrange_.size();
        } else if (t == EntityType::Face) {
            return frange_[d]->size();
        } else {
            throw std::runtime_error(
                "StructuredBase::size: Unknown entity type t");
        }
        return 0;
    }

    /// @brief Returns the total number of elements in the mesh for a
    ///        Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    size_t size(const EntityType t, const Dir d) const
    {
        return size(t, static_cast<size_t>(d));
    }

    /// @brief Returns the number of elements in the mesh for a
    ///        Cubism::EntityType for each dimension
    ///
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    MultiIndex getSize(const EntityType t, const size_t d = 0) const
    {
        if (t == EntityType::Cell) {
            return crange_.getExtent();
        } else if (t == EntityType::Node) {
            return nrange_.getExtent();
        } else if (t == EntityType::Face) {
            return frange_[d]->getExtent();
        } else {
            throw std::runtime_error(
                "StructuredBase::getSize: Unknown entity type t");
        }
        return MultiIndex();
    }

    /// @brief Returns the number of elements in the mesh for a
    ///        Cubism::EntityType for each dimension
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    MultiIndex getSize(const EntityType t, const Dir d) const
    {
        return getSize(t, static_cast<size_t>(d));
    }

    /// @brief Returns the index range for a Cubism::EntityType
    ///
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    IndexRangeType getIndexRange(const EntityType t, const size_t d = 0) const
    {
        if (t == EntityType::Cell) {
            return crange_;
        } else if (t == EntityType::Node) {
            return nrange_;
        } else if (t == EntityType::Face) {
            return *frange_[d];
        } else {
            throw std::runtime_error(
                "StructuredBase::getIndexRange: Unknown entity type t");
        }
        return IndexRangeType();
    }

    /// @brief Returns the index range for a Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    IndexRangeType getIndexRange(const EntityType t, const Dir d) const
    {
        return getIndexRange(t, static_cast<size_t>(d));
    }

    /// @brief Returns the physical mesh extent in each dimension
    PointType getExtent() const { return range_.getExtent(); }
    /// @brief Returns the volume of the mesh
    RealType getVolume() const { return range_.getVolume(); }
    /// @brief Returns the local origin of the mesh
    PointType getOrigin() const { return range_.getBegin(); }
    /// @brief Returns the global origin of the mesh
    PointType getGlobalOrigin() const { return global_origin_; }
    /// @brief Returns the physical domain range spanned by the mesh
    RangeType getRange() const { return range_; }

    /// @brief Returns true if this is a sub-mesh (subset of some main mesh)
    bool isSubMesh() const { return type_ == MeshHull::SubMesh; }

    /// @brief Converts a local flat index to a local multi-dimensional index
    ///        for a Cubism::EntityType
    ///
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    MultiIndex
    getMultiIndex(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getMultiIndex_(i, t, d);
    }

    /// @brief Converts a local flat index to a local multi-dimensional index
    ///        for a Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    MultiIndex
    getMultiIndex(const size_t i, const EntityType t, const Dir d) const
    {
        return getMultiIndex_(i, t, static_cast<size_t>(d));
    }

    /// @brief Converts a local flat index to a global multi-dimensional index
    ///        for a Cubism::EntityType
    ///
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    MultiIndex
    getGlobalIndex(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getGlobalIndex(getMultiIndex_(i, t, d), t, d);
    }

    /// @brief Converts a local flat index to a global multi-dimensional index
    ///        for a Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    MultiIndex
    getGlobalIndex(const size_t i, const EntityType t, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getGlobalIndex(getMultiIndex_(i, t, di), t, di);
    }

    /// @brief Converts a local multi-dimensional index to a global
    ///        multi-dimensional index for a Cubism::EntityType
    ///
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    MultiIndex getGlobalIndex(const MultiIndex &p,
                              const EntityType t,
                              const size_t d = 0) const
    {
        if (t == EntityType::Cell) {
            return crange_.getBegin() + p;
        } else if (t == EntityType::Node) {
            return nrange_.getBegin() + p;
        } else if (t == EntityType::Face) {
            return frange_[d]->getBegin() + p;
        } else {
            throw std::runtime_error(
                "StructuredBase::getGlobalIndex: Unknown entity type t");
        }
        return p;
    }

    /// @brief Converts a local multi-dimensional index to a global
    ///        multi-dimensional index for a Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    MultiIndex
    getGlobalIndex(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return getGlobalIndex(p, t, static_cast<size_t>(d));
    }

    /// @brief Converts a local flat index to global coordinates for a
    ///        Cubism::EntityType
    ///
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    PointType getGlobalCoords(const size_t i,
                              const EntityType t,
                              const size_t d = 0) const
    {
        return global_origin_ + getCoords(i, t, d);
    }

    /// @brief Converts a local flat index to global coordinates for a
    ///        Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    PointType
    getGlobalCoords(const size_t i, const EntityType t, const Dir d) const
    {
        return global_origin_ + getCoords(i, t, d);
    }

    /// @brief Converts a local multi-dimensional index to global coordinates
    ///        for a Cubism::EntityType
    ///
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    PointType getGlobalCoords(const MultiIndex &p,
                              const EntityType t,
                              const size_t d = 0) const
    {
        return global_origin_ + getCoords(p, t, d);
    }

    /// @brief Converts a local multi-dimensional index to global coordinates
    ///        for a Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    PointType
    getGlobalCoords(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return global_origin_ + getCoords(p, t, d);
    }

    /// @brief Converts an entity iterator to global coordinates
    ///
    /// @param it Entity iterator
    PointType getGlobalCoords(const iterator &it) const
    {
        return global_origin_ + getCoords(it);
    }

    /// @brief Converts a local flat index to local coordinates for a
    ///        Cubism::EntityType
    ///
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    PointType
    getCoords(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getCoords_(getMultiIndex_(i, t, d), t, d);
    }

    /// @brief Converts a local flat index to local coordinates for a
    ///        Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    PointType getCoords(const size_t i, const EntityType t, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getCoords_(getMultiIndex_(i, t, di), t, di);
    }

    /// @brief Converts a local multi-dimensional index to local coordinates for
    ///        a Cubism::EntityType
    ///
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    PointType
    getCoords(const MultiIndex &p, const EntityType t, const size_t d = 0) const
    {
        return getCoords_(p, t, d);
    }

    /// @brief Converts a local multi-dimensional index to local coordinates for
    ///        a Cubism::EntityType
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param d Direction indicator (for Face types only)
    template <typename Dir>
    PointType
    getCoords(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return getCoords_(p, t, static_cast<size_t>(d));
    }

    /// @brief Converts an entity iterator to local coordinates
    ///
    /// @param it Entity iterator
    PointType getCoords(const iterator &it) const
    {
        return getCoords_(*it, it.getEntity(), it.getDirection());
    }

    /// @brief Returns cell volume for local flat index
    ///
    /// @param i Local flat index
    RealType getCellVolume(const size_t i) const
    {
        return getCellVolume_(getMultiIndex_(i, EntityType::Cell));
    }

    /// @brief Returns cell volume for local multi-dimensional index
    ///
    /// @param p Local multi-dimensional index
    RealType getCellVolume(const MultiIndex &p) const
    {
        return getCellVolume_(p);
    }

    /// @brief Returns cell volume for an entity iterator
    ///
    /// @param it Entity iterator
    RealType getCellVolume(const iterator &it) const
    {
        return getCellVolume_(*it);
    }

    /// @brief Returns cell size along all dimensions for local flat index
    ///
    /// @param i Local flat index
    PointType getCellSize(const size_t i) const
    {
        return getCellSize_(getMultiIndex_(i, EntityType::Cell));
    }

    /// @brief Returns cell size along all dimensions for local
    ///        multi-dimensional index
    ///
    /// @param p Local multi-dimensional index
    PointType getCellSize(const MultiIndex &p) const { return getCellSize_(p); }

    /// @brief Returns cell size along all dimensions for an entity iterator
    ///
    /// @param it Entity iterator
    PointType getCellSize(const iterator &it) const
    {
        return getCellSize_(*it);
    }

    /// @brief Returns the surface vector for a face.  Vector points outwards of
    ///        specified reference cell
    ///
    /// @param i Local flat index of face
    /// @param ci Local flat index of reference cell (determines orientation)
    /// @param d Direction of face
    PointType getSurface(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d);
    }

    /// @brief Returns the surface vector for a face.  Vector points outwards of
    ///        specified reference cell
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index of face
    /// @param ci Local flat index of reference cell (determines orientation)
    /// @param d Direction of face
    template <typename Dir>
    PointType getSurface(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di);
    }

    /// @brief Returns the surface vector for a face.  Vector points outwards of
    ///        specified reference cell
    ///
    /// @param p Local multi-dimensional index of face
    /// @param ci Local multi-dimensional index of reference cell (determines
    ///           orientation)
    /// @param d Direction of face
    PointType
    getSurface(const MultiIndex &p, const MultiIndex &ci, const size_t d) const
    {
        return getSurface_(p, ci, d);
    }

    /// @brief Returns the surface vector for a face.  Vector points outwards of
    ///        specified reference cell
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param p Local multi-dimensional index of face
    /// @param ci Local multi-dimensional index of reference cell (determines
    ///           orientation)
    /// @param d Direction of face
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

    /// @brief Returns the surface area for a face relative to a reference cell
    ///
    /// @param i Local flat index of face
    /// @param ci Local flat index of reference cell
    /// @param d Direction of face
    RealType
    getSurfaceArea(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d)
            .norm();
    }

    /// @brief Returns the surface area for a face relative to a reference cell
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index of face
    /// @param ci Local flat index of reference cell
    /// @param d Direction of face
    template <typename Dir>
    RealType getSurfaceArea(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di)
            .norm();
    }

    /// @brief Returns the surface area for a face relative to a reference cell
    ///
    /// @param p Local multi-dimensional index of face
    /// @param ci Local multi-dimensional index of reference cell
    /// @param d Direction of face
    RealType getSurfaceArea(const MultiIndex &p,
                            const MultiIndex &ci,
                            const size_t d) const
    {
        return getSurface_(p, ci, d).norm();
    }

    /// @brief Returns the surface area for a face relative to a reference cell
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param p Local multi-dimensional index of face
    /// @param ci Local multi-dimensional index of reference cell
    /// @param d Direction of face
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

    /// @brief Returns the surface normal vector for a face.  Vector points
    ///        outwards of specified reference cell
    ///
    /// @param i Local flat index of face
    /// @param ci Local flat index of reference cell (determines orientation)
    /// @param d Direction of face
    PointType
    getSurfaceNormal(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d)
            .unit();
    }

    /// @brief Returns the surface normal vector for a face.  Vector points
    ///        outwards of specified reference cell
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param i Local flat index of face
    /// @param ci Local flat index of reference cell (determines orientation)
    /// @param d Direction of face
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

    /// @brief Returns the surface normal vector for a face.  Vector points
    ///        outwards of specified reference cell
    ///
    /// @param p Local multi-dimensional index of face
    /// @param ci Local multi-dimensional index of reference cell (determines
    ///           orientation)
    /// @param d Direction of face
    PointType getSurfaceNormal(const MultiIndex &p,
                               const MultiIndex &ci,
                               const size_t d) const
    {
        return getSurface_(p, ci, d).unit();
    }

    /// @brief Returns the surface normal vector for a face.  Vector points
    ///        outwards of specified reference cell
    ///
    /// @tparam Dir Special type that defines a cast to size_t
    /// @param p Local multi-dimensional index of face
    /// @param ci Local multi-dimensional index of reference cell (determines
    ///           orientation)
    /// @param d Direction of face
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
    const MeshHull type_;   // mesh hull type (full mesh or sub-mesh)
    const MeshClass class_; // class category of mesh
    const RangeType range_; // range of mesh domain in physical space
    const PointType global_origin_; // origin in global mesh
    const IndexRangeType crange_;   // index range spanned by cells
    const IndexRangeType nrange_;   // index range spanned by nodes
    std::vector<const IndexRangeType *>
        frange_; // index ranges spanned by faces

    /// @brief Pure virtual method to obtain local coordinates
    ///
    /// @param p Local multi-dimensional index
    /// @param t Cubism::EntityType
    /// @param dir Entity direction (for faces only)
    virtual PointType getCoords_(const MultiIndex &p,
                                 const EntityType t,
                                 const size_t dir) const = 0;

    /// @brief Pure virtual method to obtain cell volume
    ///
    /// @param p Local multi-dimensional index
    virtual RealType getCellVolume_(const MultiIndex &p) const = 0;

    /// @brief Pure virtual method to obtain cell size
    ///
    /// @param p Local multi-dimensional index
    virtual PointType getCellSize_(const MultiIndex &p) const = 0;

    /// @brief Pure virtual method to obtain surface vector
    ///
    /// @param p Local multi-dimensional index
    /// @param ci Local multi-dimensional index of reference cell (determines
    ///           orientation)
    /// @param dir Entity direction (for faces only)
    virtual PointType getSurface_(const MultiIndex &p,
                                  const MultiIndex &ci,
                                  const size_t dir) const = 0;

private:
    void initFaceRange_(const IndexRangeType &r)
    {
        for (size_t i = 0; i < frange_.size(); ++i) {
            frange_[i] = new IndexRangeType(
                r.getBegin(), r.getEnd() + MultiIndex::getUnitVector(i));
        }
    }

    void copyFaceRange_(const StructuredBase &m)
    {
        for (size_t i = 0; i < frange_.size(); ++i) {
            frange_[i] = new IndexRangeType(*m.frange_[i]);
        }
    }

    void disposeFaceRange_()
    {
        for (auto fr : frange_) {
            if (fr) {
                delete fr;
            }
        }
    }

    MultiIndex
    getMultiIndex_(const size_t i, const EntityType t, const size_t d = 0) const
    {
        if (t == EntityType::Cell) {
            return crange_.getMultiIndex(i);
        } else if (t == EntityType::Node) {
            return nrange_.getMultiIndex(i);
        } else if (t == EntityType::Face) {
            return frange_[d]->getMultiIndex(i);
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
