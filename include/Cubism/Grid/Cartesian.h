// File       : Cartesian.h
// Created    : Sun Jan 05 2020 04:34:54 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian grid composed of block fields
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIAN_H_QBSFTWK7
#define CARTESIAN_H_QBSFTWK7

#include "Alloc/AlignedBlockAllocator.h"
#include "Block/Field.h"
#include "Core/Range.h"

#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)

template <typename TData,
          typename TMesh,
          typename TEntity = Cubism::EntityType::Cell,
          size_t RANK = 0,
          template <typename> class TAlloc = AlignedBlockAllocator>
class Cartesian
{
public:
    using MeshType = TMesh;
    using IndexRangeType = typename MeshType::IndexRangeType;
    using MultiIndex = typename IndexRangeType::MultiIndex;
    using RangeType = typename MeshType::RangeType;
    using PointType = typename MeshType::PointType;

    struct FieldState {
        size_t rank;
        size_t comp;
        MultiIndex idx;
        MeshType *mesh;
    };

private:
    template <typename TE>
    struct ScalarFieldBase {
        using Type =
            Block::Field<Block::Data<TData, TE, MeshType::Dim>, FieldState>;
    };

    template <>
    struct ScalarFieldBase<Cubism::EntityType::Face> {
        using Type = Block::FaceFieldAll<TData, FieldState, MeshType::Dim>;
    };

    template <size_t R>
    struct TensorFieldBase {
        using Type =
            Block::TensorField<typename ScalarFieldBase<TEntity>::Type, R>;
    };

    template <>
    struct TensorFieldBase<0> {
        using Type = typename ScalarFieldBase<TEntity>::Type;
    };

protected:
    using MeshHull = typename MeshType::MeshHull;

public:
    using FieldType = TensorFieldBase<RANK>;
    using FieldView = Block::FieldView<FieldType>;
    using FieldContainer = Block::FieldContainer<FieldView>;

    Cartesian(const MultiIndex &nblocks,
              const MultiIndex &block_cells,
              const PointType &start = PointType(0),
              const PointType &end = PointType(1),
              const PointType &gorigin = PointType(0))
        : nblocks_(nblocks), block_cells_(block_cells)
    {
        alloc_();
        global_mesh_ = new MeshType(gorigin,
                                    RangeType(start, end),
                                    IndexRangeType(nblocks * block_cells),
                                    MeshHull::FullMesh);

        const PointType block_extent =
            global_mesh_.getExtent() / PointType(nblocks_);
        const IndexRangeType block_range(nblocks_);

        char *src = static_cast<char *>(data_);
        for (size_t i = 0; i < block_range.size(); ++i) {
            // initialize the field state
            field_states_.push_back(FieldState{0});

            // compute block mesh
            const MultiIndex bi = block_range.getMultiIndex(i);
            const PointType bstart =
                global_mesh_.getOrigin() + PointType(bi) * block_extent;
            const PointType bend = bstart + block_extent;
            const MIndex cells = block_cells_;
            MIndex nodes = cells;
            std::vector<IndexRangeType> face_ranges;
            for (size_t d = 0; d < MeshType::Dim; ++d) {
                MIndex faces(cells);
                if (bi[d] == nblocks_[d] - 1) {
                    ++nodes[i];
                    ++faces[i];
                }
                face_ranges.push_back(IndexRangeType(faces));
            }
            const IndexRangeType cell_range(cells);
            const IndexRangeType node_range(nodes);
            field_meshes_.push_back(new MeshType(gorigin,
                                                 RangeType(bstart, bend),
                                                 cell_range,
                                                 node_range,
                                                 face_ranges,
                                                 MeshHull::SubMesh));

            std::vector<IndexRangeType> rl;
            std::vector<DataType *> dl;
            std::vector<size_t> bl;
            std::vector<FieldState *> sl;
            for (size_t c = 0; c < FieldType::NComponents; ++c) {
                sl.push_back(&field_states_.back());
                dl.push_back(static_cast<DataType *>(
                    src + c * component_bytes_ + i * block_bytes_));
                bl.push_back(block_bytes_);
                if (BlockData::EntityType == EntityType::Cell) {
                    rl.push_back(cell_range);
                } else if (BlockData::EntityType == EntityType::Node) {
                    rl.push_back(node_range);
                } else if (BlockData::EntityType == EntityType::Face) {
                    // XXX: [fabianw@mavt.ethz.ch; 2020-01-05] does not work
                    // like that!
                    rl.push_back(face_ranges[c]);
                }
            }

            if (FieldType::NComponents > 1) {
                fields_.pushBack(FieldView(rl, dl, bl, sl));
            } else {
                fields_.pushBack(FieldView(rl[0], dl[0], bl[0], sl[0]));
            }
        }
    }
    virtual ~Cartesian() { dispose_(); }

protected:
    MeshType *global_mesh_;
    FieldContainer fields_;
    std::vector<FieldState> field_states_;
    std::vector<MeshType *> field_meshes_;

private:
    using BlockData = typename FieldType::BaseType;
    using DataType = typename FieldType::DataType;

    const MultiIndex nblocks_;
    const MultiIndex block_cells_;
    DataType *data_;
    size_t block_elements_;
    size_t nfaces_;
    size_t block_bytes_;
    size_t component_bytes_;
    size_t all_bytes_;
    TAlloc<DataType> blk_alloc_;

    void alloc_()
    {
        if (BlockData::EntityType == EntityType::Cell) {
            block_elements_ = block_cells_.prod();
        } else if (BlockData::EntityType == EntityType::Node) {
            block_elements_ = (block_cells_ + MultiIndex(1)).prod();
        } else if (BlockData::EntityType == EntityType::Face) {
            // XXX: [fabianw@mavt.ethz.ch; 2020-01-05] slightly more than needed
            // but easier for alignment
            block_elements_ = (block_cells_ + MultiIndex(1)).prod();
        }
        all_bytes_ = block_elements_ * sizeof(DataType);

        // align at CUBISM_ALIGNMENT byte boundary
        all_bytes_ = ((all_bytes_ + CUBISM_ALIGNMENT - 1) / CUBISM_ALIGNMENT) *
                     CUBISM_ALIGNMENT;

        // aligned block bytes (may be larger than the minimum number of bytes
        // needed)
        block_bytes_ = all_bytes_;

        // number of bytes for a single component slice of all blocks in the
        // Cartesian topology
        component_bytes_ = block_bytes_ * nblocks_.prod();

        // if the EntityType of this Cartesian grid is Face, we need to take
        // that into account for the allocated data
        nfaces_ = 1;
        if (BlockData::EntityType == EntityType::Face) {
            nfaces_ = MeshType::Dim;
        }

        // total number of bytes
        all_bytes_ = nfaces_ * component_bytes_ * FieldType::NComponents;

        // get the allocation
        assert(all_bytes_ > 0);
        blk_alloc_.allocate(all_bytes_);
    }

    void dealloc_() { blk_alloc_.deallocate(data_); }

    void dispose_()
    {
        if (global_mesh_) {
            delete global_mesh_;
        }
        for (auto fm : field_meshes_) {
            if (fm) {
                delete fm;
            }
        }
        dealloc_();
    }
};

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* CARTESIAN_H_QBSFTWK7 */
