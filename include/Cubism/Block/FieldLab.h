// File       : FieldLab.h
// Created    : Mon Feb 10 2020 06:53:39 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Data laboratory with stencil specification
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef FIELDLAB_H_8PPKFFY4
#define FIELDLAB_H_8PPKFFY4

#include "Cubism/Block/Data.h"
#include "Cubism/Block/FieldLabLoader.h"
#include "Cubism/Common.h"
#include <cassert>
#include <cstring>
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

/**
 * @brief Field laboratory
 * @tparam TField Field type to map to the lab
 *
 * @rst
 * A ``FieldLab`` is an extended data structure to include ghost cells for a
 * given stencil. Loading a lab takes care of loading the ghost cells from
 * neighboring block fields and applies boundary conditions if present.  The
 * default is periodic if no boundary conditions are specified otherwise.
 * @endrst
 * */
template <typename TField>
class FieldLab
    // Alternatively inheritance could be from TField::BlockDataType which would
    // not allow to use a different data allocator.  The inheritance below
    // allows to change Cubism::AlignedBlockAllocator for a FieldLab if needed
    // at some point.
    : public Data<typename TField::DataType,
                  TField::EntityType,
                  TField::IndexRangeType::Dim,
                  Cubism::AlignedBlockAllocator<typename TField::DataType>>
{
    using BaseType =
        Data<typename TField::DataType,
             TField::EntityType,
             TField::IndexRangeType::Dim,
             Cubism::AlignedBlockAllocator<typename TField::DataType>>;
    using BCVector = typename TField::BCVector;
    using LabLoader =
        Block::FieldLabLoader<TField, BaseType::IndexRangeType::Dim>;
    using BoolVec = typename LabLoader::BoolVec;
    using STDFunction = typename LabLoader::STDFunction;

    using BaseType::blk_alloc_;
    using BaseType::block_;
    using BaseType::bytes_;
    using BaseType::range_;

public:
    using BlockDataType = BaseType;
    using typename BaseType::BlockBytes;
    using typename BaseType::DataType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using Index = typename MultiIndex::DataType;
    using StencilType = typename LabLoader::StencilType;
    using FieldType = TField;

    /** @brief Main constructor */
    FieldLab()
        : BaseType(IndexRangeType()), is_allocated_(false),
          block_data_(nullptr), field_(nullptr), lab_begin_(0)
    {
    }

    FieldLab(const FieldLab &c) = delete;
    FieldLab(FieldLab &&c) = delete;
    FieldLab &operator=(const FieldLab &c) = delete;
    FieldLab &operator=(FieldLab &&c) = delete;
    ~FieldLab() override = default;

    using iterator = Core::MultiIndexIterator<IndexRangeType::Dim>;
    iterator begin() noexcept { return iterator(loader_.curr_range, 0); }
    iterator begin() const noexcept { return iterator(loader_.curr_range, 0); }
    iterator end() noexcept
    {
        return iterator(loader_.curr_range, loader_.curr_range.size());
    }
    iterator end() const noexcept
    {
        return iterator(loader_.curr_range, loader_.curr_range.size());
    }

    /**
     * @brief Allocate data lab memory block for a given stencil
     * @param s Target stencil
     * @param max_request_range Maximum index range to be processed in the lab
     * @param force Force a reallocation for subsequent requests
     *
     * @rst
     * The method may be called subsequently.with a different stencil and
     * maximum range.  If the new request can be processed using the current
     * allocation no new memory is allocated unless the ``force`` flag is true.
     * @endrst
     */
    void allocate(const StencilType &s,
                  const IndexRangeType &max_request_range,
                  const bool force = false)
    {
        // TODO: [fabianw@mavt.ethz.ch; 2020-02-11] Should offsets be const? any
        // runtime difference?

        // 1. Assign new stencil
        // 2. Compute full lab extent (incl. halos)
        // 3. Clear existing allocation and allocate aligned lab block

        // 1.
        loader_.curr_stencil = s;

        // 2.
        // add two extra memory locations in each direction for equal treatment
        // of block fields in a grid topology which may differ by one cell for
        // boundary adjacent blocks (depends on Cubism::EntityType).
        MultiIndex max_extent = max_request_range.getExtent() + 2;
        // maximum range that can be handled by this lab allocation
        max_range_ = IndexRangeType(max_extent);

        // memory alignment adjustment
        const typename MultiIndex::DataType n_per_align =
            CUBISM_ALIGNMENT / sizeof(DataType);

        // number of extra computational elements needed due to stencil
        lab_begin_ = -loader_.curr_stencil.getBegin();
        const MultiIndex lab_end = loader_.curr_stencil.getEnd() - 1;

        // adjust fastest moving index to alignment requirement
        lab_begin_[0] =
            ((lab_begin_[0] + n_per_align - 1) / n_per_align) * n_per_align;

        // expand the maximum range extent by additional ghost cells on the far
        // end and adjust the fastest moving index to alignment requirement
        max_extent += lab_end;
        max_extent[0] =
            ((max_extent[0] + n_per_align - 1) / n_per_align) * n_per_align;

        // total lab extent including ghost cells and alignment requirements
        const MultiIndex lab_extent = lab_begin_ + max_extent;

        // if this is a subsequent call, can we recycle the current allocation?
        const bool can_reuse = lab_extent <= range_.getExtent();

        // this is the new lab extent for this allocation call
        range_ = IndexRangeType(lab_extent);

        // if there is an existing allocation and no forced re-allocation
        // needed, exit here
        if (!force && is_allocated_ && can_reuse) {
            block_data_ = block_ + range_.getFlatIndex(lab_begin_);
            return;
        }

        // 3.
        BaseType::deallocBlock_();
        BaseType::allocBlock_();
        block_data_ = block_ + range_.getFlatIndex(lab_begin_);
        is_allocated_ = true;
    }

    /**
     * @brief Lab data loader
     * @param fid Multi-dimensional index of target block field
     * @param id2field Index mapping function for block fields
     * @param apply_bc Flag whether to apply boundary conditions
     * @param extern_bc Pointer to external boundary conditions
     *
     * @rst
     * The ``id2field`` mapping function takes a multi-dimensional block field
     * index as an argument and returns a reference to the corresponding block
     * field.  The function must map indices periodically.
     * @endrst
     */
    template <typename Functor = STDFunction>
    void loadData(const MultiIndex &fid,
                  Functor &id2field,
                  const bool apply_bc = true,
                  const BCVector *extern_bc = nullptr)
    {
        static_assert(FieldType::Class == Cubism::FieldClass::Scalar,
                      "FieldLab: field class must be scalar.");
        if (!is_allocated_) {
            throw std::runtime_error(
                "FieldLab: can not load lab data when not allocated first");
        }

        // 1. load the block field data
        // 2. load the halos
        // 3. apply boundary conditions

        // 1.
        field_ = &(id2field(fid));
        loader_.curr_range = field_->getIndexRange();
        loader_.curr_labrange = IndexRangeType(
            loader_.curr_stencil.getBegin(),
            loader_.curr_range.getExtent() + loader_.curr_stencil.getEnd() - 1);
        loader_.loadInner(*field_, block_, range_, lab_begin_);

        // 2.
        const BCVector &bcs = (extern_bc) ? *extern_bc : field_->getBC();
        BoolVec periodic(true);
        MultiIndex skip(1);
        for (const auto bc : bcs) {
            const auto info = bc->getBoundaryInfo();
            assert(info.dir < IndexRangeType::Dim);
            periodic[info.dir] = info.is_periodic;
            skip[info.dir] = (info.side == 0) ? -1 : 1;
        }
        loader_.loadGhosts(
            fid, id2field, block_, range_, lab_begin_, periodic, skip);

        // 3.
        if (apply_bc) {
            for (const auto bc : bcs) {
                (*bc)(*this);
            }
        }
    }

    /**
     * @brief Lab data loader
     * @param fid Multi-dimensional index of target block field
     * @param id2field Index mapping function for block fields
     * @param boundaries Vector of boundary conditions
     * @param apply_bc Flag whether to apply boundary conditions
     *
     * @rst
     * This loader is a convenience wrapper around the default loader if
     * boundary conditions need to be enforced.  This wrapper applies the
     * boundary conditions specified in ``boundaries`` instead of the boundary
     * conditions specified in the block field with index ``fid`` (if any).  The
     * vector ``boundaries`` may be empty.
     * @endrst
     */
    template <typename Functor = STDFunction>
    void loadData(const MultiIndex &fid,
                  Functor &id2field,
                  const BCVector &boundaries,
                  const bool apply_bc = true)
    {
        loadData(fid, id2field, apply_bc, &boundaries);
    }

    /**
     * @brief Linear data access
     * @param p Local multi-dimensional index
     * @return Reference to data element
     *
     * @rst
     * The local index ``p`` may reference halo cells.  For example, ``p{-1,0}``
     * indexes the first halo cell in the ``x`` direction that is adjacent to
     * the inner domain.
     * @endrst
     */
    DataType &operator[](const MultiIndex &p)
    {
        assert(range_.isIndex(p + lab_begin_));
        if (1 == IndexRangeType::Dim) {
            return block_[p[0] + lab_begin_[0]];
        } else if (2 == IndexRangeType::Dim) {
            return block_[p[0] + lab_begin_[0] +
                          range_.sizeDim(0) * (p[1] + lab_begin_[1])];
        } else if (3 == IndexRangeType::Dim) {
            return block_[p[0] + lab_begin_[0] +
                          range_.sizeDim(0) *
                              (p[1] + lab_begin_[1] +
                               range_.sizeDim(1) * (p[2] + lab_begin_[2]))];
        } else {
            return BaseType::operator[](range_.getFlatIndex(p + lab_begin_));
        }
    }

    /**
     * @brief Linear data access
     * @param p Local multi-dimensional index
     * @return ``const`` reference to data element
     *
     * @rst
     * The local index ``p`` may reference halo cells.  For example, ``p{-1,0}``
     * indexes the first halo cell in the ``x`` direction that is adjacent to
     * the inner domain.
     * @endrst
     */
    const DataType &operator[](const MultiIndex &p) const
    {
        assert(range_.isIndex(p + lab_begin_));
        if (1 == IndexRangeType::Dim) {
            return block_[p[0] + lab_begin_[0]];
        } else if (2 == IndexRangeType::Dim) {
            return block_[p[0] + lab_begin_[0] +
                          range_.sizeDim(0) * (p[1] + lab_begin_[1])];
        } else if (3 == IndexRangeType::Dim) {
            return block_[p[0] + lab_begin_[0] +
                          range_.sizeDim(0) *
                              (p[1] + lab_begin_[1] +
                               range_.sizeDim(1) * (p[2] + lab_begin_[2]))];
        } else {
            return BaseType::operator[](range_.getFlatIndex(p + lab_begin_));
        }
    }

    /**
     * @brief Classic data access
     * @param ix Index for first dimension
     * @param iy Index for second dimension
     * @param iz Index for third dimension
     * @return Reference to data element
     *
     * This operator is only supported for dimensions 1, 2 and 3.  Accessing
     * data with this operator has less latency than the multi-index sibling.
     */
    DataType &operator()(const Index ix, const Index iy = 0, const Index iz = 0)
    {
        assert((ix >= loader_.curr_labrange.getBegin()[0]) &&
               (ix <= loader_.curr_labrange.getEnd()[0]));
        if (1 == IndexRangeType::Dim) {
            return block_[ix + lab_begin_[0]];
        } else if (2 == IndexRangeType::Dim) {
            assert((iy >= loader_.curr_labrange.getBegin()[1]) &&
                   (iy <= loader_.curr_labrange.getEnd()[1]));
            return block_[ix + lab_begin_[0] +
                          range_.sizeDim(0) * (iy + lab_begin_[1])];
        } else if (3 == IndexRangeType::Dim) {
            assert((iy >= loader_.curr_labrange.getBegin()[1]) &&
                   (iy <= loader_.curr_labrange.getEnd()[1]));
            assert((iz >= loader_.curr_labrange.getBegin()[2]) &&
                   (iz <= loader_.curr_labrange.getEnd()[2]));
            return block_[ix + lab_begin_[0] +
                          range_.sizeDim(0) *
                              (iy + lab_begin_[1] +
                               range_.sizeDim(1) * (iz + lab_begin_[2]))];
        }
        throw std::runtime_error("FieldLab: operator() not supported");
    }

    /**
     * @brief Classic data access
     * @param ix Index for first dimension
     * @param iy Index for second dimension
     * @param iz Index for third dimension
     * @return ``const`` reference to data element
     *
     * This operator is only supported for dimensions 1, 2 and 3.  Accessing
     * data with this operator has less latency than the multi-index sibling.
     */
    const DataType &
    operator()(const Index ix, const Index iy = 0, const Index iz = 0) const
    {
        assert((ix >= loader_.curr_labrange.getBegin()[0]) &&
               (ix <= loader_.curr_labrange.getEnd()[0]));
        if (1 == IndexRangeType::Dim) {
            return block_[ix + lab_begin_[0]];
        } else if (2 == IndexRangeType::Dim) {
            assert((iy >= loader_.curr_labrange.getBegin()[1]) &&
                   (iy <= loader_.curr_labrange.getEnd()[1]));
            return block_[ix + lab_begin_[0] +
                          range_.sizeDim(0) * (iy + lab_begin_[1])];
        } else if (3 == IndexRangeType::Dim) {
            assert((iy >= loader_.curr_labrange.getBegin()[1]) &&
                   (iy <= loader_.curr_labrange.getEnd()[1]));
            assert((iz >= loader_.curr_labrange.getBegin()[2]) &&
                   (iz <= loader_.curr_labrange.getEnd()[2]));
            return block_[ix + lab_begin_[0] +
                          range_.sizeDim(0) *
                              (iy + lab_begin_[1] +
                               range_.sizeDim(1) * (iz + lab_begin_[2]))];
        }
        throw std::runtime_error("FieldLab: operator() not supported");
    }

    /**
     * @brief Get pointer to inner block data
     * @return Pointer to first data element of base field
     *
     * @rst
     * The returned pointer points to the first element defined in
     * ``loader_.curr_range``.
     * @endrst
     */
    DataType *getInnerData() { return block_data_; }

    /**
     * @brief Get pointer to inner block data
     * @return Pointer to first data element of base field
     *
     * @rst
     * The returned pointer points to the first element defined in
     * ``loader_.curr_range``.
     * @endrst
     */
    const DataType *getInnerData() const { return block_data_; }

    /**
     * @brief Get currently active (loaded) stencil
     * @return ``const`` reference to ``StencilType``
     */
    const StencilType &getActiveStencil() const { return loader_.curr_stencil; }

    /**
     * @brief Get currently active (loaded) index range
     * @return ``const`` reference to ``IndexRangeType``
     */
    const IndexRangeType &getActiveRange() const { return loader_.curr_range; }

    /**
     * @brief Get reference to currently active (loaded) field
     * @return Reference to ``FieldType``
     */
    FieldType &getActiveField()
    {
        if (!field_) {
            throw std::runtime_error("FieldLab: no field loaded");
        }
        return *field_;
    }

    /**
     * @brief Get reference to currently active (loaded) field
     * @return ``const`` reference to ``FieldType``
     */
    const FieldType &getActiveField() const
    {
        if (!field_) {
            throw std::runtime_error("FieldLab: no field loaded");
        }
        return *field_;
    }

    /**
     * @brief Get currently active (loaded) lab index range
     * @return Index range including ghost indices
     *
     * @rst
     * The index range begin is identical to the active stencil begin, which is
     * ``< 0``.
     * @endrst
     */
    IndexRangeType getActiveLabRange() const { return loader_.curr_labrange; }

    /**
     * @brief Get the maximum range that the lab can hold
     * @return Index range of maximum span
     */
    IndexRangeType getMaximumRange() const { return max_range_; }

    /**
     * @brief Check if lab is allocated
     * @return True if laboratory is allocated
     */
    bool isAllocated() const { return is_allocated_; }

    /**
     * @brief Get byte utilization of block
     * @return Structure of byte usage for this instance
     */
    BlockBytes getMemoryFootprint() const override
    {
        BlockBytes bb = {};
        if (is_allocated_) {
            bb.allocated = this->bytes_;
            bb.used = (loader_.curr_range.getExtent() +
                       loader_.curr_stencil.getEnd() -
                       loader_.curr_stencil.getBegin() - 1)
                          .prod() *
                      sizeof(DataType);
        }
        return bb;
    }

protected:
    bool allocBlock_() override
    {
        // allocate nothing during constructor call
        return true;
    }

private:
    bool is_allocated_;
    IndexRangeType max_range_;
    LabLoader loader_;
    DataType *block_data_; // start of block data
    FieldType *field_;     // currently loaded field

    // offset helper
    MultiIndex lab_begin_;
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* FIELDLAB_H_8PPKFFY4 */
