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

#include <cassert>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity = Cubism::EntityType::Cell,
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
    template <Cubism::EntityType TE>
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
    /// @brief Block (tensor) field type
    using FieldType = TensorFieldBase<RANK>;
    /// @brief View type for block fields
    using FieldView = Block::FieldView<FieldType>;
    /// @brief Container type for field views
    using FieldContainer = Block::FieldContainer<FieldView>;

    static constexpr size_t Rank = RANK;
    static constexpr size_t NComponents = FieldType::NComponents;
    static constexpr typename Cubism::EntityType EntityType = TEntity;

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
        initBlockFields_(global_mesh_->getGlobalOrigin());
        initFieldViews_();

        assert(fields_.size() == field_states_.size());
        assert(fields_.size() == field_meshes_.size());
        assert(fields_.size() == tensor_fields_.size());
    }

    virtual ~Cartesian() { dispose_(); }

protected:
    /// @brief Type of components in the tensor fields
    using FieldBaseType = typename FieldType::FieldType;

    MeshType *global_mesh_;
    FieldContainer fields_;
    std::vector<FieldState> field_states_;
    std::vector<MeshType *> field_meshes_;
    std::vector<FieldType *> tensor_fields_;

private:
    using BlockData = typename FieldType::BaseType;
    using DataType = typename FieldType::DataType;

    const MultiIndex nblocks_;
    const MultiIndex block_cells_;
    DataType *data_;
    size_t block_elements_;
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
        size_t nfaces_ = 1;
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
        for (auto tf : tensor_fields_) {
            if (tf) {
                if (NComponents == 1) {
                    delete tf;
                } else {
                    for (auto c : *tf) {
                        if (c) {
                            delete c;
                        }
                    }
                    delete tf;
                }
            }
        }
        for (auto fv : fields_) {
            if (fv) {
                delete fv;
            }
        }
        dealloc_();
    }

    void initBlockFields_(const PointType &gorigin)
    {
        const PointType block_extent =
            global_mesh_->getExtent() / PointType(nblocks_);
        const IndexRangeType block_range(nblocks_);

        char *const src = static_cast<char *>(data_);
        for (size_t i = 0; i < block_range.size(); ++i) {
            // initialize the field state
            field_states_.push_back(FieldState{0});
            FieldState &fs = &field_states_.back();

            // compute block mesh
            const MultiIndex bi = block_range.getMultiIndex(i);
            const PointType bstart =
                global_mesh_->getOrigin() + PointType(bi) * block_extent;
            const PointType bend = bstart + block_extent;
            const MultiIndex cells = block_cells_;
            MultiIndex nodes = cells;
            std::vector<IndexRangeType> face_ranges;
            for (size_t d = 0; d < MeshType::Dim; ++d) {
                MultiIndex faces(cells);
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
            fs.idx = bi;
            fs.mesh = field_meshes_.back();

            // generate views
            if (NComponents == 1) { // scalar field: FieldType == FieldBaseType
                if (EntityType == Cubism::EntityType::Face) {
                    std::vector<DataType *> dl;
                    std::vector<size_t> bl;
                    std::vector<FieldState *> sl;
                    for (size_t d = 0; d < MeshType::Dim; ++d) {
                        char *dst =
                            src + i * block_bytes_ + d * component_bytes_;
                        dl.push_back(static_cast<DataType *>(dst));
                        bl.push_back(block_bytes_);
                        sl.push_back(&fs); // all point to the same state
                    }
                    FieldType *sf = new FieldType(face_ranges, dl, bl, sl);
                    tensor_fields_.push_back(sf);
                } else if (EntityType == Cubism::EntityType::Node) {
                    char *dst = src + i * block_bytes_;
                    FieldType *sf = new FieldType(node_range,
                                                  static_cast<DataType *>(dst),
                                                  block_bytes_,
                                                  &fs);
                    tensor_fields_.push_back(sf);
                } else if (EntityType == Cubism::EntityType::Cell) {
                    char *dst = src + i * block_bytes_;
                    FieldType *sf = new FieldType(cell_range,
                                                  static_cast<DataType *>(dst),
                                                  block_bytes_,
                                                  &fs);
                    tensor_fields_.push_back(sf);
                }
            } else { // tensor field: FieldType != FieldBaseType
                FieldType *tf = new FieldType(); // empty tensor field
                tensor_fields_.push_back(tf);
                for (size_t c = 0; c < NComponents; ++c) {
                    if (EntityType == Cubism::EntityType::Face) {
                        std::vector<DataType *> dl;
                        std::vector<size_t> bl;
                        std::vector<FieldState *> sl;
                        for (size_t d = 0; d < MeshType::Dim; ++d) {
                            char *dst = src +
                                        c * MeshType::Dim * component_bytes_ +
                                        d * component_bytes_ + i * block_bytes_;
                            dl.push_back(static_cast<DataType *>(dst));
                            bl.push_back(block_bytes_);
                            sl.push_back(&fs); // all point to the same state
                        }
                        FieldBaseType *sf =
                            new FieldBaseType(face_ranges, dl, bl, sl);
                        tf->pushBack(sf);
                    } else if (EntityType == Cubism::EntityType::Node) {
                        char *dst =
                            src + c * component_bytes_ + i * block_bytes_;
                        FieldBaseType *sf =
                            new FieldBaseType(node_range,
                                              static_cast<DataType *>(dst),
                                              block_bytes_,
                                              &fs);
                        tf->pushBack(sf);
                    } else if (EntityType == Cubism::EntityType::Cell) {
                        char *dst =
                            src + c * component_bytes_ + i * block_bytes_;
                        FieldBaseType *sf =
                            new FieldBaseType(cell_range,
                                              static_cast<DataType *>(dst),
                                              block_bytes_,
                                              &fs);
                        tf->pushBack(sf);
                    }
                }
            }
        }
    }

    void initFieldViews_()
    {
        for (auto tf : tensor_fields_) {
            fields_.pushBack(new FieldView(*tf));
        }
    }
};

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity,
          size_t RANK,
          template <typename>
          class TAlloc>
constexpr size_t Cartesian<TData, TMesh, TEntity, RANK, TAlloc>::Rank;

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity,
          size_t RANK,
          template <typename>
          class TAlloc>
constexpr size_t Cartesian<TData, TMesh, TEntity, RANK, TAlloc>::NComponents;

template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity,
          size_t RANK,
          template <typename>
          class TAlloc>
constexpr typename Cubism::EntityType
    Cartesian<TData, TMesh, TEntity, RANK, TAlloc>::EntityType;

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* CARTESIAN_H_QBSFTWK7 */
