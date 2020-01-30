// File       : BlockFieldAssembler.h
// Created    : Wed Jan 08 2020 12:02:57 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian block field assembler and specializations
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef BLOCKFIELDASSEMBLER_H_GXHBUCPV
#define BLOCKFIELDASSEMBLER_H_GXHBUCPV

#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"
#include <cstddef>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)

/**
 * @brief Block field assembler for an externally allocated region of memory
 * @tparam T Field data type
 * @tparam RANK Tensor rank
 * @tparam Entity Entity type
 * @tparam State Field state type
 * @tparam Mesh Mesh type of global and block meshes
 * */
template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          typename State,
          typename Mesh>
struct BlockFieldAssembler {
    /** @brief Main field type */
    using BaseType = typename Block::
        FieldTypeFactory<T, RANK, Entity, Mesh::Dim, State>::Type;
    /** @brief Container of tensor fields */
    using FieldContainer = Block::FieldContainer<BaseType>;
    /** @brief Scalar field type of components in main field */
    using FieldType = typename BaseType::FieldType;
    /** @brief State (meta) information for block */
    using FieldState = typename FieldType::FieldStateType;
    /** @brief Underlying block data manager used by the assembled fields */
    using BlockData = typename FieldType::BlockDataType;
    /** @brief Synonym type for T */
    using DataType = typename BlockData::DataType;
    /** @brief Type of mesh */
    using MeshType = Mesh;
    /** @brief Type of mesh hull (full mesh or sub-mesh) */
    using MeshIntegrity = typename MeshType::MeshIntegrity;
    /** @brief Index range type */
    using IndexRangeType = typename MeshType::IndexRangeType;
    /** @brief Type for higher dimensional index */
    using MultiIndex = typename IndexRangeType::MultiIndex;
    /** @brief Type for physical domain ranges spanned by ``MeshType`` */
    using RangeType = typename MeshType::RangeType;
    /** @brief Type of point in physical domain */
    using PointType = typename MeshType::PointType;

    /**
     * @brief State (meta data) of the assembled fields
     */
    std::vector<FieldState *> field_states;
    /**
     * @brief Block meshes corresponding to the assembled fields
     */
    std::vector<MeshType *> field_meshes;
    /**
     * @brief Assembled fields into the external memory
     */
    FieldContainer fields;

    /**
     * @brief Main assembly routine
     * @param src Pointer to the beginning of the externally allocated memory
     * @param mesh Global mesh in which contains the assembled block meshes
     * @param block_range Range of blocks for which to assemble the fields
     * @param block_cells Number of cells for individual blocks
     * @param scale Scaling factor for total number of blocks
     * @param block_bytes Number of bytes occupied by each block
     * @param component_bytes Number of bytes per tensor component (must be
     *        larger or equal to ``nblocks.prod() * block_bytes``)
     *
     * @rst
     * Assembles block fields and its sub mesh on a Cartesian topology using
     * the external data ``src``.
     * @endrst
     */
    void assemble(DataType *src,
                  const MeshType &mesh,
                  const IndexRangeType &block_range,
                  const MultiIndex &block_cells,
                  const MultiIndex &scale,
                  const size_t block_bytes,
                  const size_t component_bytes)
    {
        dispose();

        const MultiIndex nblocks = block_range.getExtent();
        const MultiIndex all_blocks = scale * nblocks;
        const PointType block_extent = mesh.getExtent() / PointType(nblocks);
        const MultiIndex c0 = mesh.getIndexRange(EntityType::Cell).getBegin();

        char *const base = reinterpret_cast<char *>(src);
        std::vector<IndexRangeType> face_ranges(MeshType::Dim);
        std::vector<IndexRangeType> A;
        std::vector<DataType *> B;
        std::vector<size_t> C;
        std::vector<FieldState *> D;
        std::vector<std::vector<IndexRangeType>> AA;
        std::vector<std::vector<DataType *>> BB;
        std::vector<std::vector<size_t>> CC;
        std::vector<std::vector<FieldState *>> DD;

        for (size_t i = 0; i < block_range.size(); ++i) {
            // initialize the field state
            FieldState *fs = new FieldState();
            field_states.push_back(fs);

            // compute block mesh
            const MultiIndex bi = block_range.getMultiIndex(i); // local index
            const MultiIndex gbi = block_range.getBegin() + bi; // global index
            const MultiIndex cstart = c0 + bi * block_cells;
            const PointType bstart =
                mesh.getOrigin() + PointType(bi) * block_extent;
            const PointType bend = bstart + block_extent;
            const MultiIndex cells = block_cells;
            MultiIndex nodes = cells;
            for (size_t d = 0; d < MeshType::Dim; ++d) {
                MultiIndex faces(cells);
                if (gbi[d] == all_blocks[d] - 1) {
                    ++nodes[d];
                    ++faces[d];
                }
                face_ranges[d] = IndexRangeType(cstart, cstart + faces);
            }
            const IndexRangeType cell_range(cstart, cstart + cells);
            const IndexRangeType node_range(cstart, cstart + nodes);
            MeshType *fm = new MeshType(mesh.getGlobalOrigin(),
                                        RangeType(bstart, bend),
                                        cell_range,
                                        node_range,
                                        face_ranges,
                                        MeshIntegrity::SubMesh);
            assert(fm != nullptr);
            field_meshes.push_back(fm);
            fs->idx = bi;
            fs->mesh = fm;

            // generate fields
            A.clear();
            B.clear();
            C.clear();
            D.clear();
            AA.clear();
            BB.clear();
            CC.clear();
            DD.clear();
            if (Entity == Cubism::EntityType::Face) {
                for (size_t d = 0; d < MeshType::Dim; ++d) {
                    for (size_t c = 0; c < BaseType::NComponents; ++c) {
                        char *dst =
                            base + d * BaseType::NComponents * component_bytes +
                            c * component_bytes + i * block_bytes;
                        A.push_back(face_ranges[d]);
                        B.push_back(reinterpret_cast<DataType *>(dst));
                        C.push_back(block_bytes);
                        D.push_back(fs);
                    }
                    AA.push_back(A);
                    BB.push_back(B);
                    CC.push_back(C);
                    DD.push_back(D);
                    A.clear();
                    B.clear();
                    C.clear();
                    D.clear();
                }
            } else if (Entity == Cubism::EntityType::Node) {
                for (size_t c = 0; c < BaseType::NComponents; ++c) {
                    char *dst = base + c * component_bytes + i * block_bytes;
                    A.push_back(node_range);
                    B.push_back(reinterpret_cast<DataType *>(dst));
                    C.push_back(block_bytes);
                    D.push_back(fs);
                }
                AA.push_back(A);
                BB.push_back(B);
                CC.push_back(C);
                DD.push_back(D);
            } else if (Entity == Cubism::EntityType::Cell) {
                for (size_t c = 0; c < BaseType::NComponents; ++c) {
                    char *dst = base + c * component_bytes + i * block_bytes;
                    A.push_back(cell_range);
                    B.push_back(reinterpret_cast<DataType *>(dst));
                    C.push_back(block_bytes);
                    D.push_back(fs);
                }
                AA.push_back(A);
                BB.push_back(B);
                CC.push_back(C);
                DD.push_back(D);
            } else {
                throw std::runtime_error(
                    "BlockFieldAssembler: Unknown entity type");
            }
            fields.pushBack(new BaseType(AA, BB, CC, DD));
        }
    }

    /**
     * @brief Dispose assembled fields
     */
    void dispose()
    {
        fields.clear();
        for (auto fm : field_meshes) {
            if (fm) {
                delete fm;
            }
        }
        for (auto fs : field_states) {
            if (fs) {
                delete fs;
            }
        }
        field_meshes.clear();
        field_states.clear();
    }
};

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* BLOCKFIELDASSEMBLER_H_GXHBUCPV */
