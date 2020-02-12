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
#include <cstring>
#include <functional>
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

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

    using BaseType::blk_alloc_;
    using BaseType::block_;
    using BaseType::bytes_;
    using BaseType::range_;

public:
    using typename BaseType::DataType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using StencilType = Core::Stencil<IndexRangeType::Dim>;

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

    void allocate(const StencilType &s,
                  const IndexRangeType &max_range,
                  const bool force = false)
    {
        // TODO: [fabianw@mavt.ethz.ch; 2020-02-11] does fastest moving index
        // need to be pitched for SIMD registers?
        //
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

    void loadData(const MultiIndex &fid,
                  const ID2Field &id2field,
                  const double time = 0,
                  const bool apply_bc = true)
    {
        static_assert(TField::Class == Cubism::FieldClass::Scalar,
                      "DataLab: field class must be scalar.");
        if (!is_allocated_) {
            throw std::runtime_error(
                "DataLab: can not load lab data when not allocated first.");
        }

        // 1. load the block field data
        // 2. load the halos / apply boundary conditions

        // 1.
        const TField &f0 = id2field(fid);
        active_range_ = f0.getIndexRange();

        // TODO: [fabianw@mavt.ethz.ch; 2020-02-12] specialized versions
        // DataType *dst = block_data_;
        // const DataType *src = f0.getData();
        // const size_t l0 = range_.sizeDim(0);        // lab stride
        // const size_t n0 = active_range_.sizeDim(0); // data stride
        // const size_t bytes0 = n0 * sizeof(DataType);
        // const size_t N = active_range_.size();

        const auto range_end = active_range_.end();
        for (auto it = active_range_.begin(); it != range_end; ++it) {
            *(block_ + range_.getFlatIndex(*it + lab_begin_)) =
                f0[it.getFlatIndex()];
        }

        // 2.
        auto boundaries = f0.getBC(); // get list of boundary conditions

        // XXX: [fabianw@mavt.ethz.ch; 2020-02-12] continue periodic BC

        //{
        //    const bool xperiodic = is_xperiodic();
        //    const bool yperiodic = is_yperiodic();
        //    const bool zperiodic = is_zperiodic();

        //    const bool xskin = info.index[0]==0 ||
        //    info.index[0]==grid.getBlocksPerDimension(0)-1; const bool yskin =
        //    info.index[1]==0 ||
        //    info.index[1]==grid.getBlocksPerDimension(1)-1; const bool zskin =
        //    info.index[2]==0 ||
        //    info.index[2]==grid.getBlocksPerDimension(2)-1;

        //    const int xskip = info.index[0]==0 ? -1 : 1;
        //    const int yskip = info.index[1]==0 ? -1 : 1;
        //    const int zskip = info.index[2]==0 ? -1 : 1;

        //    for(int icode=0; icode<27; icode++)
        //    {
        //        if (icode == 1*1 + 3*1 + 9*1) continue;

        //        const int code[3] = { icode%3-1, (icode/3)%3-1,
        //        (icode/9)%3-1};

        //        if (!xperiodic && code[0] == xskip && xskin) continue;
        //        if (!yperiodic && code[1] == yskip && yskin) continue;
        //        if (!zperiodic && code[2] == zskip && zskin) continue;

        //        if (!istensorial && abs(code[0])+abs(code[1])+abs(code[2])>1)
        //        continue;

        //        const int s[3] = {
        //            code[0]<1? (code[0]<0 ? m_stencilStart[0]:0 ) : nX,
        //            code[1]<1? (code[1]<0 ? m_stencilStart[1]:0 ) : nY,
        //            code[2]<1? (code[2]<0 ? m_stencilStart[2]:0 ) : nZ };

        //        const int e[3] = {
        //            code[0]<1? (code[0]<0 ? 0:nX ) : nX+m_stencilEnd[0]-1,
        //            code[1]<1? (code[1]<0 ? 0:nY ) : nY+m_stencilEnd[1]-1,
        //            code[2]<1? (code[2]<0 ? 0:nZ ) : nZ+m_stencilEnd[2]-1};

        //        if (!grid.avail(info.index[0] + code[0], info.index[1] +
        //        code[1], info.index[2] + code[2])) continue;

        //        BlockType& b = grid(info.index[0] + code[0], info.index[1] +
        //        code[1], info.index[2] + code[2]);

        //#if 1
        //        const int m_vSize0 = m_cacheBlock->getSize(0); //m_vSize[0];
        //        const int m_nElemsPerSlice =
        //        m_cacheBlock->getNumberOfElementsPerSlice();
        //        //m_nElementsPerSlice;

        //        const int my_ix = s[0]-m_stencilStart[0];

        //        //printf("iy : %d %d\n", s[1], e[1]);
        //        const int bytes = (e[0]-s[0])*sizeof(ElementType);
        //        if (!bytes) continue;
        //        for(int iz=s[2]; iz<e[2]; iz++)
        //        {
        //            const int my_izx = (iz-m_stencilStart[2])*m_nElemsPerSlice
        //            + my_ix; #if 0 for(int iy=s[1]; iy<e[1]; iy++)
        //            {
        //                #if 1	// ...
        //                //char * ptrDest =
        //                (char*)&m_cacheBlock->Access(s[0]-m_stencilStart[0],
        //                iy-m_stencilStart[1], iz-m_stencilStart[2]); char *
        //                ptrDest = (char*)&m_cacheBlock->LinAccess(my_izx +
        //                (iy-m_stencilStart[1])*m_vSize0);

        //                const char * ptrSrc = (const char*)&b(s[0] -
        //                code[0]*BlockType::sizeX, iy -
        //                code[1]*BlockType::sizeY, iz -
        //                code[2]*BlockType::sizeZ); memcpy2((char *)ptrDest,
        //                (char *)ptrSrc, bytes); #else for(int ix=s[0];
        //                ix<e[0]; ix++)
        //                    m_cacheBlock->Access(ix-m_stencilStart[0],
        //                    iy-m_stencilStart[1], iz-m_stencilStart[2]) =
        //                    (ElementType)b(ix - code[0]*BlockType::sizeX, iy -
        //                    code[1]*BlockType::sizeY, iz -
        //                    code[2]*BlockType::sizeZ);
        //                #endif
        //            }
        //            #else
        //            if ((e[1]-s[1]) % 4 != 0)
        //            {
        //                for(int iy=s[1]; iy<e[1]; iy++)
        //                {
        //                    char * ptrDest =
        //                    (char*)&m_cacheBlock->LinAccess(my_izx +
        //                    (iy-m_stencilStart[1])*m_vSize0);

        //                    const char * ptrSrc = (const char*)&b(s[0] -
        //                    code[0]*BlockType::sizeX, iy -
        //                    code[1]*BlockType::sizeY, iz -
        //                    code[2]*BlockType::sizeZ); const int cpybytes =
        //                    (e[0]-s[0])*sizeof(ElementType); memcpy2((char
        //                    *)ptrDest, (char *)ptrSrc, cpybytes);
        //                }
        //            }
        //            else
        //            {
        //                for(int iy=s[1]; iy<e[1]; iy+=4)
        //                {
        //                    char * ptrDest0 =
        //                    (char*)&m_cacheBlock->LinAccess(my_izx +
        //                    (iy+0-m_stencilStart[1])*m_vSize0); char *
        //                    ptrDest1 = (char*)&m_cacheBlock->LinAccess(my_izx
        //                    + (iy+1-m_stencilStart[1])*m_vSize0); char *
        //                    ptrDest2 = (char*)&m_cacheBlock->LinAccess(my_izx
        //                    + (iy+2-m_stencilStart[1])*m_vSize0); char *
        //                    ptrDest3 = (char*)&m_cacheBlock->LinAccess(my_izx
        //                    + (iy+3-m_stencilStart[1])*m_vSize0);

        //                    const char * ptrSrc0 = (const char*)&b(s[0] -
        //                    code[0]*BlockType::sizeX, iy + 0 -
        //                    code[1]*BlockType::sizeY, iz -
        //                    code[2]*BlockType::sizeZ); const char * ptrSrc1 =
        //                    (const char*)&b(s[0] - code[0]*BlockType::sizeX,
        //                    iy + 1 - code[1]*BlockType::sizeY, iz -
        //                    code[2]*BlockType::sizeZ); const char * ptrSrc2 =
        //                    (const char*)&b(s[0] - code[0]*BlockType::sizeX,
        //                    iy + 2 - code[1]*BlockType::sizeY, iz -
        //                    code[2]*BlockType::sizeZ); const char * ptrSrc3 =
        //                    (const char*)&b(s[0] - code[0]*BlockType::sizeX,
        //                    iy + 3 - code[1]*BlockType::sizeY, iz -
        //                    code[2]*BlockType::sizeZ);

        //                    memcpy2((char *)ptrDest0, (char *)ptrSrc0, bytes);
        //                    memcpy2((char *)ptrDest1, (char *)ptrSrc1, bytes);
        //                    memcpy2((char *)ptrDest2, (char *)ptrSrc2, bytes);
        //                    memcpy2((char *)ptrDest3, (char *)ptrSrc3, bytes);
        //                }
        //            }
        //            #endif
        //        }
        //#else
        //        const int off_x = - code[0]*nX + m_stencilStart[0];

        if (apply_bc) {
            for (auto bc : boundaries) {
                (*bc)(*this, time);
            }
        }
    }

    // template <typename TField>
    // void loadData(const MultiIndex &fid,
    //               const ID2Field<TField> &id2field,
    //               const std::vector<Cubism::BC::Base*> &boundary_conditions,
    //               const double time = 0,
    //               const bool apply_bc = true)

protected:
    bool allocBlock_() override { return true; }

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
