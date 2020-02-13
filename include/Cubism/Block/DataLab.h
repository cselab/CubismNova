// File       : DataLab.h
// Created    : Mon Feb 10 2020 06:53:39 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Data laboratory with stencil specification
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef DATALAB_H_O091Y6A2
#define DATALAB_H_O091Y6A2

#include "Cubism/Block/Data.h"
#include "Cubism/Common.h"
#include "Cubism/Core/Stencil.h"
#include "Cubism/Core/Vector.h"
#include "Cubism/Math.h"
#include <cassert>
#include <cstring>
#include <functional>
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

/**
 * @brief Data laboratory
 * @tparam TField Field type to map to the lab
 *
 * @rst
 * A ``DataLab`` is an extended data structure to include ghost cells for a
 * given stencil. Loading a lab takes care of loading the ghost cells from
 * neighboring block fields and applies boundary conditions if present.  The
 * default is periodic if no boundary conditions are specified otherwise.
 * @endrst
 * */
template <typename TField>
class DataLab
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

    using ID2Field =
        std::function<const TField &(const typename BaseType::MultiIndex &)>;
    using BoolVec = Core::Vector<bool, TField::IndexRangeType::Dim>;

    using BaseType::blk_alloc_;
    using BaseType::block_;
    using BaseType::bytes_;
    using BaseType::range_;

public:
    using typename BaseType::BlockBytes;
    using typename BaseType::DataType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using StencilType = Core::Stencil<IndexRangeType::Dim>;

    /** @brief Main constructor */
    explicit DataLab()
        : BaseType(IndexRangeType()), is_allocated_(false),
          active_stencil_(0, 1), active_range_(IndexRangeType()),
          block_data_(nullptr), lab_begin_(0), lab_end_(0)
    {
    }

    DataLab(const DataLab &c) = delete;
    DataLab(DataLab &&c) = delete;
    DataLab &operator=(const DataLab &c) = delete;
    DataLab &operator=(DataLab &&c) = delete;
    ~DataLab() override = default;

    using iterator = Core::MultiIndexIterator<IndexRangeType::Dim>;
    iterator begin() noexcept { return iterator(active_range_, 0); }
    iterator begin() const noexcept { return iterator(active_range_, 0); }
    iterator end() noexcept
    {
        return iterator(active_range_, active_range_.size());
    }
    iterator end() const noexcept
    {
        return iterator(active_range_, active_range_.size());
    }

    /**
     * @brief Allocate data lab memory block for a given stencil
     * @param s Target stencil
     * @param max_range Maximum index range to be processed in the lab
     * @param force Force a reallocation for subsequent requests
     *
     * @rst
     * The method may be called subsequently.with a different stencil and
     * maximum range.  If the new request can be processed using the current
     * allocation no new memory is allocated unless the ``force`` flag is true.
     * @endrst
     */
    void allocate(const StencilType &s,
                  const IndexRangeType &max_range,
                  const bool force = false)
    {
        // TODO: [fabianw@mavt.ethz.ch; 2020-02-11] Should offsets be const? any
        // runtime difference?

        // 1. Assign new stencil
        // 2. Compute full lab extent (incl. halos)
        // 3. Clear existing allocation and allocate aligned lab block

        // 1.
        active_stencil_ = s;

        // 2.
        // add two extra memory locations in each direction for equal treatment
        // of block fields in a grid topology which may differ by one cell for
        // boundary adjacent blocks (depends on Cubism::EntityType).
        MultiIndex max_extent = max_range.getExtent() + 2;

        const typename MultiIndex::DataType n_per_align =
            CUBISM_ALIGNMENT / sizeof(DataType);
        lab_begin_ = -active_stencil_.getBegin();
        lab_end_ = active_stencil_.getEnd() - 1;
        lab_begin_[0] =
            ((lab_begin_[0] + n_per_align - 1) / n_per_align) * n_per_align;
        max_extent += lab_end_;
        max_extent[0] =
            ((max_extent[0] + n_per_align - 1) / n_per_align) * n_per_align;
        const MultiIndex lab_extent = lab_begin_ + max_extent;
        const bool can_reuse = lab_extent <= range_.getExtent();
        range_ = IndexRangeType(lab_extent);
        if (!force && is_allocated_ && can_reuse) {
            block_data_ = block_ + range_.getFlatIndex(lab_begin_);
            return;
        }

        // 3.
        BaseType::deallocBlock_();
        if (!BaseType::allocBlock_()) {
            throw std::runtime_error(
                "DataLab: error allocating new lab block memory.");
        }
        block_data_ = block_ + range_.getFlatIndex(lab_begin_);
        is_allocated_ = true;
    }

    /**
     * @brief Lab data loader
     * @param fid Multi-dimensional index of target block field
     * @param id2field Index mapping function for block fields
     * @param apply_bc Flag whether to apply boundary conditions
     *
     * @rst
     * The ``id2field`` mapping function takes a multi-dimensional block field
     * index as an argument and returns a reference to the corresponding block
     * field.  The function must map indices periodically.
     * @endrst
     */
    void loadData(const MultiIndex &fid,
                  const ID2Field &id2field,
                  const bool apply_bc = true)
    {
        static_assert(TField::Class == Cubism::FieldClass::Scalar,
                      "DataLab: field class must be scalar.");
        if (!is_allocated_) {
            throw std::runtime_error(
                "DataLab: can not load lab data when not allocated first.");
        }

        // 1. load the block field data
        // 2. load the halos
        // 3. apply boundary conditions

        // 1.
        const TField &f0 = id2field(fid);
        active_range_ = f0.getIndexRange();

        const auto range_end = active_range_.end();
        for (auto it = active_range_.begin(); it != range_end; ++it) {
            *(block_ + range_.getFlatIndex(*it + lab_begin_)) =
                f0[it.getFlatIndex()];
        }

        // 2.
        auto boundaries = f0.getBC(); // get list of boundary conditions
        BoolVec periodic(true);
        MultiIndex skip(1);
        for (auto bc : boundaries) {
            const auto info = bc->getBoundaryInfo();
            assert(info.dir < IndexRangeType::Dim);
            periodic[info.dir] = info.is_periodic;
            skip[info.dir] = (info.side == 0) ? -1 : 1;
        }

        const MultiIndex one(1);
        const IndexRangeType nbr_range(0, 3);
        const size_t neighbors = nbr_range.size();
        const size_t me = neighbors / 2;
        const MultiIndex active_extent = active_range_.getExtent();
        const MultiIndex halo_extent =
            active_extent + active_stencil_.getEnd() - one;
        const MultiIndex stencil_begin = active_stencil_.getBegin();
        for (size_t i = 0; i < neighbors; ++i) {
            if (i == me) {
                continue;
            }
            const MultiIndex bi = nbr_range.getMultiIndex(i) - one;

            typename MultiIndex::DataType isum = 0;
            bool skip_current = false;
            for (size_t j = 0; j < IndexRangeType::Dim; ++j) {
                if (!periodic[j] && bi[j] == skip[j]) {
                    skip_current = true;
                    break;
                }
                isum += Cubism::myAbs(bi[j]);
            }
            if (skip_current) {
                continue;
            }
            if (!active_stencil_.isTensorial() && isum > 1) {
                continue;
            }

            MultiIndex begin;
            MultiIndex end;
            for (size_t j = 0; j < IndexRangeType::Dim; ++j) {
                begin[j] = (bi[j] < 1) ? ((bi[j] < 0) ? stencil_begin[j] : 0)
                                       : active_extent[j];
                end[j] = (bi[j] < 1) ? ((bi[j] < 0) ? 0 : active_extent[j])
                                     : halo_extent[j];
            }
            const IndexRangeType halo_range(begin, end);
            const MultiIndex lab_begin = begin + lab_begin_;
            const MultiIndex nbr_begin = begin - bi * active_extent;

            const auto &f1 = id2field(fid + bi);
            for (auto &p : halo_range) {
                *(block_ + range_.getFlatIndex(p + lab_begin)) =
                    f1[p + nbr_begin];
            }
        }

        // 3.
        if (apply_bc) {
            for (auto bc : boundaries) {
                (*bc)(*this);
            }
        }
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
        return BaseType::operator[](range_.getFlatIndex(p + lab_begin_));
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
        return BaseType::operator[](range_.getFlatIndex(p + lab_begin_));
    }

    /**
     * @brief Get pointer to inner block data
     * @return Pointer to first data element of base field
     *
     * @rst
     * The returned pointer points to the first element defined in
     * ``active_range_``.
     * @endrst
     */
    DataType *getInnerData() { return block_data_; }

    /**
     * @brief Get pointer to inner block data
     * @return Pointer to first data element of base field
     *
     * @rst
     * The returned pointer points to the first element defined in
     * ``active_range_``.
     * @endrst
     */
    const DataType *getInnerData() const { return block_data_; }

    /**
     * @brief Get currently active stencil
     * @return ``const`` reference to ``StencilType``
     */
    const StencilType &getActiveStencil() const { return active_stencil_; }

    /**
     * @brief Get currently active index range
     * @return ``const`` reference to ``IndexRangeType``
     */
    const IndexRangeType &getActiveRange() const { return active_range_; }

    /**
     * @brief Get currently active lab index range
     * @return Index range including ghost indices
     *
     * @rst
     * The index range begin is identical to the active stencil begin, which is
     * ``< 0``.
     * @endrst
     */
    IndexRangeType getActiveLabRange() const
    {
        return IndexRangeType(active_stencil_.getBegin(),
                              active_range_.getExtent() +
                                  active_stencil_.getEnd() - MultiIndex(1));
    }

    /**
     * @brief Get byte utilization of block
     * @return Structure of byte usage for this instance
     */
    BlockBytes getMemoryFootprint() const override
    {
        BlockBytes bb = {};
        if (is_allocated_) {
            bb.allocated = this->bytes_;
            bb.used = (active_range_.getExtent() + active_stencil_.getEnd() -
                       active_stencil_.getBegin() - MultiIndex(1))
                          .prod() *
                      sizeof(DataType);
        }
        return bb;
    }

    // template <typename TField>
    // void loadData(const MultiIndex &fid,
    //               const ID2Field<TField> &id2field,
    //               const std::vector<Cubism::BC::Base*> &boundary_conditions,
    //               const double time = 0,
    //               const bool apply_bc = true)

    // TODO: [fabianw@mavt.ethz.ch; 2020-02-12] specialized versions (notes)
    // DataType *dst = block_data_;
    // const DataType *src = f0.getData();
    // const size_t l0 = range_.sizeDim(0);        // lab stride
    // const size_t n0 = active_range_.sizeDim(0); // data stride
    // const size_t bytes0 = n0 * sizeof(DataType);
    // const size_t N = active_range_.size();

protected:
    bool allocBlock_() override
    {
        // allocate nothing during constructor call
        return true;
    }

private:
    bool is_allocated_;
    StencilType active_stencil_;
    IndexRangeType active_range_;
    DataType *block_data_; // start of block data

    // lab memory
    MultiIndex lab_begin_;
    MultiIndex lab_end_;
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* DATALAB_H_O091Y6A2 */
