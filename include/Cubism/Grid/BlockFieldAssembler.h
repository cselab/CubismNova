// File       : BlockFieldAssembler.h
// Created    : Wed Jan 08 2020 12:02:57 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian block field assembler and specializations
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef BLOCKFIELDASSEMBLER_H_GXHBUCPV
#define BLOCKFIELDASSEMBLER_H_GXHBUCPV

#include "Block/Field.h"
#include "Common.h"
#include <cstddef>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)

/// @brief Scalar field type factory
///
/// @tparam TEntity Cubism::EntityType
/// @tparam TData Field data type
/// @tparam TState Field state type
/// @tparam DIM Field dimension
template <Cubism::EntityType TEntity,
          typename TData,
          typename TState,
          size_t DIM>
struct ScalarFieldBase {
    using Type =
        Cubism::Block::Field<Cubism::Block::Data<TData, TEntity, DIM>, TState>;
};

/// @brief Scalar field type factory (specialization for FaceFieldAll types)
///
/// @tparam TData Field data type
/// @tparam TState Field state type
/// @tparam DIM Field dimension
template <typename TData, typename TState, size_t DIM>
struct ScalarFieldBase<Cubism::EntityType::Face, TData, TState, DIM> {
    using Type = Cubism::Block::FaceFieldAll<TData, TState, DIM>;
};

/// @brief Tensor field type factory
///
/// @tparam RANK Tensor rank
/// @tparam TEntity Cubism::EntityType
/// @tparam TData Field data type
/// @tparam TState Field state type
/// @tparam DIM Field dimension
template <size_t RANK,
          Cubism::EntityType TEntity,
          typename TData,
          typename TState,
          size_t DIM>
struct TensorFieldBase {
    using Type = Cubism::Block::TensorField<
        typename ScalarFieldBase<TEntity, TData, TState, DIM>::Type,
        RANK>;
};

/// @brief Tensor field type factory (specialization for scalar fields)
///
/// @tparam TEntity Cubism::EntityType
/// @tparam TData Field data type
/// @tparam TState Field state type
/// @tparam DIM Field dimension
template <Cubism::EntityType TEntity,
          typename TData,
          typename TState,
          size_t DIM>
struct TensorFieldBase<0, TEntity, TData, TState, DIM> {
    using Type = typename ScalarFieldBase<TEntity, TData, TState, DIM>::Type;
};

/// @brief Block field assembler for an externally allocated region of memory
///
/// @tparam TEntity Cubism::EntityType
/// @tparam TFData Field data type
/// @tparam TFState Field state type
/// @tparam TMesh Mesh type of global and block meshes
/// @tparam RANK Tensor rank of field
template <Cubism::EntityType TEntity,
          typename TFData,
          typename TFState,
          typename TMesh,
          size_t RANK>
struct BlockFieldAssembler {
    /// @brief Main tensor field type
    using FieldType =
        typename TensorFieldBase<RANK, TEntity, TFData, TFState, TMesh::Dim>::
            Type;
    /// @brief Container of tensor fields
    using FieldContainer = Block::FieldContainer<FieldType>;
    /// @brief Field type of components in main tensor field
    using FieldBaseType = typename FieldType::FieldType;
    /// @brief State (meta) information for block
    using FieldState = TFState;
    /// @brief Underlying block data manager used by the assembled fields
    using BlockData = typename FieldType::BaseType;
    /// @brief Synonym type for TFData
    using DataType = typename BlockData::DataType;
    /// @brief Type of mesh
    using MeshType = TMesh;
    /// @brief Type of mesh hull (full mesh or sub-mesh)
    using MeshHull = typename MeshType::MeshHull;
    /// @brief Index range type
    using IndexRangeType = typename MeshType::IndexRangeType;
    /// @brief Type for higher dimensional index
    using MultiIndex = typename IndexRangeType::MultiIndex;
    /// @brief Type for physical domain ranges spanned by MeshType
    using RangeType = typename MeshType::RangeType;
    /// @brief Type of point in physical domain (a MeshType::Dim vector of
    ///        float or double)
    using PointType = typename MeshType::PointType;

    /// @brief State (meta data) of the assembled fields
    std::vector<FieldState> field_states;
    /// @brief Block meshes corresponding to the assembled fields
    std::vector<MeshType *> field_meshes;
    /// @brief Assembled field views into the external memory
    FieldContainer tensor_fields;

    /// @brief Main assembly routine
    ///
    /// @param src Pointer to the beginning of the externally allocated memory
    /// @param mesh Global mesh in which contains the assembled block meshes
    /// @param nblocks Number of blocks to assemble
    /// @param block_cells Number of cells for individual blocks
    /// @param block_bytes Number of bytes occupied by each block
    /// @param component_bytes Number of bytes per tensor component (must be
    ///        larger or equal to nblocks.prod() * block_bytes)
    void assemble(TFData *src,
                  const TMesh &mesh,
                  const MultiIndex &nblocks,
                  const MultiIndex &block_cells,
                  const size_t block_bytes,
                  const size_t component_bytes)
    {
        dispose();

        const PointType block_extent = mesh.getExtent() / PointType(nblocks);
        const IndexRangeType block_range(nblocks);

        char *const base = reinterpret_cast<char *>(src);
        for (size_t i = 0; i < block_range.size(); ++i) {
            // initialize the field state
            field_states.push_back(FieldState());
            FieldState &fs = field_states.back();

            // compute block mesh
            const MultiIndex bi = block_range.getMultiIndex(i); // local index
            const PointType bstart =
                mesh.getOrigin() + PointType(bi) * block_extent;
            const PointType bend = bstart + block_extent;
            const MultiIndex cells = block_cells;
            MultiIndex nodes = cells;
            std::vector<IndexRangeType> face_ranges;
            for (size_t d = 0; d < MeshType::Dim; ++d) {
                MultiIndex faces(cells);
                if (bi[d] == nblocks[d] - 1) {
                    ++nodes[i];
                    ++faces[i];
                }
                face_ranges.push_back(IndexRangeType(faces));
            }
            const IndexRangeType cell_range(cells);
            const IndexRangeType node_range(nodes);
            field_meshes.push_back(new MeshType(mesh.getGlobalOrigin(),
                                                RangeType(bstart, bend),
                                                cell_range,
                                                node_range,
                                                face_ranges,
                                                MeshHull::SubMesh));
            fs.rank = RANK;
            fs.idx = bi;
            fs.mesh = field_meshes.back();

            // generate views
            FieldType *tf = new FieldType(); // empty tensor field
            tensor_fields.pushBack(tf);
            for (size_t c = 0; c < FieldType::NComponents; ++c) {
                fs.comp = c;
                if (TEntity == Cubism::EntityType::Face) {
                    std::vector<DataType *> dl;
                    std::vector<size_t> bl;
                    std::vector<FieldState *> sl;
                    for (size_t d = 0; d < MeshType::Dim; ++d) {
                        char *dst = base + c * MeshType::Dim * component_bytes +
                                    d * component_bytes + i * block_bytes;
                        dl.push_back(reinterpret_cast<DataType *>(dst));
                        bl.push_back(block_bytes);
                        sl.push_back(&fs); // all point to the same state
                    }
                    FieldBaseType *sf =
                        new FieldBaseType(face_ranges, dl, bl, sl);
                    tf->pushBack(sf);
                } else if (TEntity == Cubism::EntityType::Node) {
                    char *dst = base + c * component_bytes + i * block_bytes;
                    FieldBaseType *sf =
                        new FieldBaseType(node_range,
                                          reinterpret_cast<DataType *>(dst),
                                          block_bytes,
                                          &fs);
                    tf->pushBack(sf);
                } else if (TEntity == Cubism::EntityType::Cell) {
                    char *dst = base + c * component_bytes + i * block_bytes;
                    FieldBaseType *sf =
                        new FieldBaseType(cell_range,
                                          reinterpret_cast<DataType *>(dst),
                                          block_bytes,
                                          &fs);
                    tf->pushBack(sf);
                }
            }
        }
    }

    /// @brief Dispose assembled fields
    void dispose()
    {
        for (auto fm : field_meshes) {
            if (fm) {
                delete fm;
            }
        }
        tensor_fields.clear();
        field_states.clear();
        field_meshes.clear();
    }
};

/// @brief Block field assembler for an externally allocated region of memory
///
/// @tparam TEntity Cubism::EntityType
/// @tparam TFData Field data type
/// @tparam TFState Field state type
/// @tparam TMesh Mesh type of global and block meshes
///
/// This is a specialization for scalar fields (RANK = 0)
template <Cubism::EntityType TEntity,
          typename TFData,
          typename TFState,
          typename TMesh>
struct BlockFieldAssembler<TEntity, TFData, TFState, TMesh, 0> {
    /// @brief Main scalar field type
    using FieldType =
        typename TensorFieldBase<0, TEntity, TFData, TFState, TMesh::Dim>::Type;
    /// @brief Container of scalar fields
    using FieldContainer = Block::FieldContainer<FieldType>;
    /// @brief Field type of components in main scalar field
    using FieldBaseType = FieldType;
    /// @brief State (meta) information for block
    using FieldState = TFState;
    /// @brief Underlying block data manager used by the assembled fields
    using BlockData = typename FieldType::BaseType;
    /// @brief Synonym type for TFData
    using DataType = typename BlockData::DataType;
    /// @brief Type of mesh
    using MeshType = TMesh;
    /// @brief Type of mesh hull (full mesh or sub-mesh)
    using MeshHull = typename MeshType::MeshHull;
    /// @brief Index range type
    using IndexRangeType = typename MeshType::IndexRangeType;
    /// @brief Type for higher dimensional index
    using MultiIndex = typename IndexRangeType::MultiIndex;
    /// @brief Type for physical domain ranges spanned by MeshType
    using RangeType = typename MeshType::RangeType;
    /// @brief Type of point in physical domain (a MeshType::Dim vector of
    ///        float or double)
    using PointType = typename MeshType::PointType;

    /// @brief State (meta data) of the assembled fields
    std::vector<FieldState> field_states;
    /// @brief Block meshes corresponding to the assembled fields
    std::vector<MeshType *> field_meshes;
    /// @brief Assembled field views into the external memory
    FieldContainer tensor_fields;

    /// @brief Main assembly routine
    ///
    /// @param src Pointer to the beginning of the externally allocated memory
    /// @param mesh Global mesh in which contains the assembled block meshes
    /// @param nblocks Number of blocks to assemble
    /// @param block_cells Number of cells for individual blocks
    /// @param block_bytes Number of bytes occupied by each block
    /// @param component_bytes Number of bytes per tensor component (must be
    ///        larger or equal to nblocks.prod() * block_bytes)
    void assemble(TFData *src,
                  const TMesh &mesh,
                  const MultiIndex &nblocks,
                  const MultiIndex &block_cells,
                  const size_t block_bytes,
                  const size_t component_bytes)
    {
        dispose();

        const PointType block_extent = mesh.getExtent() / PointType(nblocks);
        const IndexRangeType block_range(nblocks);

        char *const base = reinterpret_cast<char *>(src);
        for (size_t i = 0; i < block_range.size(); ++i) {
            // initialize the field state
            field_states.push_back(FieldState());
            FieldState &fs = field_states.back();

            // compute block mesh
            const MultiIndex bi = block_range.getMultiIndex(i); // local index
            const PointType bstart =
                mesh.getOrigin() + PointType(bi) * block_extent;
            const PointType bend = bstart + block_extent;
            const MultiIndex cells = block_cells;
            MultiIndex nodes = cells;
            std::vector<IndexRangeType> face_ranges;
            for (size_t d = 0; d < MeshType::Dim; ++d) {
                MultiIndex faces(cells);
                if (bi[d] == nblocks[d] - 1) {
                    ++nodes[i];
                    ++faces[i];
                }
                face_ranges.push_back(IndexRangeType(faces));
            }
            const IndexRangeType cell_range(cells);
            const IndexRangeType node_range(nodes);
            field_meshes.push_back(new MeshType(mesh.getGlobalOrigin(),
                                                RangeType(bstart, bend),
                                                cell_range,
                                                node_range,
                                                face_ranges,
                                                MeshHull::SubMesh));
            fs.rank = 0;
            fs.comp = 0;
            fs.idx = bi;
            fs.mesh = field_meshes.back();

            // generate views
            if (TEntity == Cubism::EntityType::Face) {
                std::vector<DataType *> dl;
                std::vector<size_t> bl;
                std::vector<FieldState *> sl;
                for (size_t d = 0; d < MeshType::Dim; ++d) {
                    char *dst = base + i * block_bytes + d * component_bytes;
                    dl.push_back(reinterpret_cast<DataType *>(dst));
                    bl.push_back(block_bytes);
                    sl.push_back(&fs); // all point to the same state
                }
                FieldType *sf = new FieldType(face_ranges, dl, bl, sl);
                tensor_fields.pushBack(sf);
            } else if (TEntity == Cubism::EntityType::Node) {
                char *dst = base + i * block_bytes;
                FieldType *sf = new FieldType(node_range,
                                              reinterpret_cast<DataType *>(dst),
                                              block_bytes,
                                              &fs);
                tensor_fields.pushBack(sf);
            } else if (TEntity == Cubism::EntityType::Cell) {
                char *dst = base + i * block_bytes;
                FieldType *sf = new FieldType(cell_range,
                                              reinterpret_cast<DataType *>(dst),
                                              block_bytes,
                                              &fs);
                tensor_fields.pushBack(sf);
            }
        }
    }

    /// @brief Dispose assembled fields
    void dispose()
    {
        for (auto fm : field_meshes) {
            if (fm) {
                delete fm;
            }
        }
        tensor_fields.clear();
        field_states.clear();
        field_meshes.clear();
    }
};

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* BLOCKFIELDASSEMBLER_H_GXHBUCPV */
