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
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename BlockAlloc = AlignedBlockAllocator<T>>
class DataLab : public Data<T, Entity, DIM, BlockAlloc>
{
    using BaseType = Data<T, Entity, DIM, BlockAlloc>;

    using BaseType::blk_alloc_;
    using BaseType::block_;
    using BaseType::bytes_;
    using BaseType::range_;

public:
    using typename BaseType::DataType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using StencilType = Core::Stencil<DIM>;

    DataLab(const IndexRangeType &max_range) : max_range_(max_range_) {}

    ~DataLab() override = default;

    void allocate(const StencilType &s)
    {
        if (is_locked_) {
            throw std::runtime_error(
                "DataLab: can not allocate new lab when locked.");
        }

        // TODO: [fabianw@mavt.ethz.ch; 2020-02-11] does fastest moving index
        // need to be pitched for SIMD registers?
        //
        // TODO: [fabianw@mavt.ethz.ch; 2020-02-11] Should offsets be const? any
        // runtime difference?

        stencil_ = s;
        this->deallocBlock_();

        const MultiIndex full_extent = max_range_.getExtent() -
                                       stencil_.getBegin() + stencil_.getEnd() -
                                       1;

        range_ = IndexRangeType(full_extent);
        this->allocBlock_();

        is_allocated_ = true;
    }

    void loadData(/* Field, functional for neighbor fields */)
    {
        // fills the memory
        is_locked_ = true;
    }

    void unlock() { is_locked_ = false; }

private:
    IndexRangeType max_range_;
    IndexRangeType active_range_;
    bool is_allocated_;
    bool is_locked_;

    StencilType stencil_;
    DataType *data_origin_; // start of block data

    size_t x_stride_pad;    // padded x-stride
    size_t slice_stride;
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* DATALAB_H_O091Y6A2 */
