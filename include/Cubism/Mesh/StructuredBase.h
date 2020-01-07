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

    /// @brief
    ///
    /// @param gorigin
    /// @param range
    /// @param crange
    /// @param type
    /// @param cl
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

    template <typename Dir>
    EntityIterator getIterator(const EntityType t, const Dir d)
    {
        return getIterator(t, static_cast<size_t>(d));
    }

    // Basic mesh interface
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

    template <typename Dir>
    size_t size(const EntityType t, const Dir d) const
    {
        return size(t, static_cast<size_t>(d));
    }

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

    template <typename Dir>
    MultiIndex getSize(const EntityType t, const Dir d) const
    {
        return getSize(t, static_cast<size_t>(d));
    }

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

    template <typename Dir>
    IndexRangeType getIndexRange(const EntityType t, const Dir d) const
    {
        return getIndexRange(t, static_cast<size_t>(d));
    }

    PointType getExtent() const { return range_.getExtent(); }
    RealType getVolume() const { return range_.getVolume(); }
    PointType getOrigin() const { return range_.getBegin(); }
    PointType getGlobalOrigin() const { return global_origin_; }
    RangeType getRange() const { return range_; }

    bool isSubMesh() const { return type_ == MeshHull::SubMesh; }

    MultiIndex
    getMultiIndex(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getMultiIndex_(i, t, d);
    }

    template <typename Dir>
    MultiIndex
    getMultiIndex(const size_t i, const EntityType t, const Dir d) const
    {
        return getMultiIndex_(i, t, static_cast<size_t>(d));
    }

    MultiIndex
    getGlobalIndex(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getGlobalIndex(getMultiIndex_(i, t, d), t, d);
    }

    template <typename Dir>
    MultiIndex
    getGlobalIndex(const size_t i, const EntityType t, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getGlobalIndex(getMultiIndex_(i, t, di), t, di);
    }

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

    template <typename Dir>
    MultiIndex
    getGlobalIndex(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return getGlobalIndex(p, t, static_cast<size_t>(d));
    }

    PointType getGlobalCoords(const size_t i,
                              const EntityType t,
                              const size_t d = 0) const
    {
        return global_origin_ + getCoords(i, t, d);
    }

    template <typename Dir>
    PointType
    getGlobalCoords(const size_t i, const EntityType t, const Dir d) const
    {
        return global_origin_ + getCoords(i, t, d);
    }

    PointType getGlobalCoords(const MultiIndex &p,
                              const EntityType t,
                              const size_t d = 0) const
    {
        return global_origin_ + getCoords(p, t, d);
    }

    template <typename Dir>
    PointType
    getGlobalCoords(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return global_origin_ + getCoords(p, t, d);
    }

    PointType getGlobalCoords(const iterator &it) const
    {
        return global_origin_ + getCoords(it);
    }

    PointType
    getCoords(const size_t i, const EntityType t, const size_t d = 0) const
    {
        return getCoords_(getMultiIndex_(i, t, d), t, d);
    }

    template <typename Dir>
    PointType getCoords(const size_t i, const EntityType t, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getCoords_(getMultiIndex_(i, t, di), t, di);
    }

    PointType
    getCoords(const MultiIndex &p, const EntityType t, const size_t d = 0) const
    {
        return getCoords_(p, t, d);
    }

    template <typename Dir>
    PointType
    getCoords(const MultiIndex &p, const EntityType t, const Dir d) const
    {
        return getCoords_(p, t, static_cast<size_t>(d));
    }

    PointType getCoords(const iterator &it) const
    {
        return getCoords_(*it, it.getEntity(), it.getDirection());
    }

    RealType getCellVolume(const size_t i) const
    {
        return getCellVolume_(getMultiIndex_(i, EntityType::Cell));
    }

    RealType getCellVolume(const MultiIndex &pi) const
    {
        return getCellVolume_(pi);
    }

    RealType getCellVolume(const iterator &it) const
    {
        return getCellVolume_(*it);
    }

    PointType getCellSize(const size_t i) const
    {
        return getCellSize_(getMultiIndex_(i, EntityType::Cell));
    }

    PointType getCellSize(const MultiIndex &p) const { return getCellSize_(p); }

    PointType getCellSize(const iterator &it) const
    {
        return getCellSize_(*it);
    }

    PointType getSurface(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d);
    }

    template <typename Dir>
    PointType getSurface(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di);
    }

    PointType
    getSurface(const MultiIndex &p, const MultiIndex &ci, const size_t d) const
    {
        return getSurface_(p, ci, d);
    }

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

    RealType
    getSurfaceArea(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d)
            .norm();
    }

    template <typename Dir>
    RealType getSurfaceArea(const size_t i, const size_t ci, const Dir d) const
    {
        const size_t di = static_cast<size_t>(d);
        return getSurface_(getMultiIndex_(i, EntityType::Face, di),
                           getMultiIndex_(ci, EntityType::Cell),
                           di)
            .norm();
    }

    RealType getSurfaceArea(const MultiIndex &p,
                            const MultiIndex &ci,
                            const size_t d) const
    {
        return getSurface_(p, ci, d).norm();
    }

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

    PointType
    getSurfaceNormal(const size_t i, const size_t ci, const size_t d) const
    {
        return getSurface_(getMultiIndex_(i, EntityType::Face, d),
                           getMultiIndex_(ci, EntityType::Cell),
                           d)
            .unit();
    }

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

    PointType getSurfaceNormal(const MultiIndex &p,
                               const MultiIndex &ci,
                               const size_t d) const
    {
        return getSurface_(p, ci, d).unit();
    }

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

    virtual PointType getCoords_(const MultiIndex &p,
                                 const EntityType t,
                                 const size_t dir) const = 0;
    virtual RealType getCellVolume_(const MultiIndex &p) const = 0;
    virtual PointType getCellSize_(const MultiIndex &p) const = 0;
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
