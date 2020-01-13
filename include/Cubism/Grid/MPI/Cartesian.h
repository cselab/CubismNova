// File       : Cartesian.h
// Created    : Sat Jan 11 2020 01:24:00 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian MPI grid composed of block fields
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIAN_H_YIX8J1RA
#define CARTESIAN_H_YIX8J1RA

#include "Core/Vector.h"
#include "Grid/Cartesian.h"
#include <cassert>
#include <mpi.h>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)
NAMESPACE_BEGIN(MPI)

/// @brief Cartesian MPI block (tensor) field
///
/// @tparam TData Field data type
/// @tparam TMesh Mesh type to be associated with fields
/// @tparam TEntity Cubism::EntityType
/// @tparam RANK Rank of (tensor) fields
/// @tparam TAlloc Allocator for field data
template <typename TData,
          typename TMesh,
          Cubism::EntityType TEntity = Cubism::EntityType::Cell,
          size_t RANK = 0,
          template <typename> class TAlloc = AlignedBlockAllocator>
class Cartesian
    : public Cubism::Grid::Cartesian<TData, TMesh, TEntity, RANK, TAlloc>
{
    using BaseType =
        Cubism::Grid::Cartesian<TData, TMesh, TEntity, RANK, TAlloc>;
    using IntVec = typename Core::Vector<int, BaseType::Dim>;

    using BaseType::block_cells_;
    using BaseType::block_range_;
    using BaseType::nblocks_;

public:
    using typename BaseType::DataType;
    using typename BaseType::FieldContainer;
    using typename BaseType::FieldState;
    using typename BaseType::FieldType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MeshType;
    using typename BaseType::MultiIndex;
    using typename BaseType::PointType;
    using typename BaseType::RangeType;
    using typename BaseType::RealType;

    /// @brief Main constructor for a Cartesian MPI block field topology
    ///
    /// @param comm MPI communicator for this Cartesian topology
    /// @param nprocs Number of MPI processes in each dimension
    /// @param nblocks Number of blocks per rank
    /// @param block_cells Number of cells in each block
    /// @param start Physical origin for the full (all ranks) Cartesian grid
    ///              (lower left)
    /// @param end Physical end for the full (all ranks) Cartesian grid (top
    ///            right)
    /// @param gorigin Physical global origin
    Cartesian(const MPI_Comm &comm,
              const MultiIndex &nprocs,
              const MultiIndex &nblocks,
              const MultiIndex &block_cells,
              const PointType &start = PointType(0),
              const PointType &end = PointType(1),
              const PointType &gorigin = PointType(0))
        : BaseType(), comm_(comm), comm_cart_(MPI_COMM_NULL), nprocs_(nprocs)
    {
        nblocks_ = nblocks;
        block_cells_ = block_cells;

        // MPI topology
        int size;
        MPI_Comm_size(comm_, &size);
        if (size != nprocs_.prod()) {
            throw std::runtime_error("Number of processes in communicator does "
                                     "not match the requested number of ranks");
        }
        IntVec pe = IntVec(nprocs_); // number of processes
        IntVec periodic(1);
        MPI_Cart_create(comm_,
                        static_cast<int>(IntVec::Dim),
                        pe.data(),
                        periodic.data(),
                        true,
                        &comm_cart_);
        MPI_Comm_rank(comm_cart_, &rank_cart_);
        IntVec pe_index; // process index in Cartesian topology
        MPI_Cart_coords(comm_cart_,
                        rank_cart_,
                        static_cast<int>(IntVec::Dim),
                        pe_index.data());
        rank_index_ = MultiIndex(pe_index);

        // block range for this rank
        const MultiIndex bstart_rank = rank_index_ * nblocks_;
        block_range_ = IndexRangeType(bstart_rank, bstart_rank + nblocks_);

        // mesh and data topology for this rank
        const PointType extent_rank = (end - start) / PointType(nprocs_);
        const PointType start_rank =
            start + PointType(rank_index_) * extent_rank; // rank domain start
        const PointType end_rank = start_rank + extent_rank; // rank domain end
        this->initTopology_(gorigin, start_rank, end_rank, nprocs_);
    }

    Cartesian() = default;
    Cartesian(const Cartesian &c) = delete;
    Cartesian(Cartesian &&c) = delete;
    Cartesian &operator=(const Cartesian &c) = default;
    Cartesian &operator=(Cartesian &&c) = delete;
    ~Cartesian() = default;

    /// @brief Returns the global number of blocks in all dimensions
    MultiIndex getGlobalSize() const { return nprocs_ * nblocks_; }
    /// @brief Returns the number of processes in the MPI topology
    MultiIndex getNumProcs() const { return nprocs_; }
    /// @brief Returns the Cartesian index for this rank
    MultiIndex getProcIndex() const { return rank_index_; }
    /// @brief Returns the rank of this process in the Cartesian communicator
    int getCartRank() const { return rank_cart_; }
    /// @brief Returns the Cartesian MPI communicator
    MPI_Comm getCartComm() const { return comm_cart_; }
    /// @brief Returns true if this rank is root
    bool isRoot() const { return (0 == rank_cart_); }

private:
    MPI_Comm comm_;      // World communicator
    MPI_Comm comm_cart_; // Cartesian communicator
    MultiIndex nprocs_;  // Number of MPI processes
    MultiIndex rank_index_; // Cartesian index of this rank
    int rank_cart_;         // Cartesian MPI rank
};

NAMESPACE_END(MPI)
NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* CARTESIAN_H_YIX8J1RA */
