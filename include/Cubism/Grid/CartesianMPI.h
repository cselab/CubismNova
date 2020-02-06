// File       : CartesianMPI.h
// Created    : Sat Jan 11 2020 01:24:00 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Cartesian MPI grid composed of block fields
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef CARTESIANMPI_H_INHL4O2K
#define CARTESIANMPI_H_INHL4O2K

#include "Cubism/Common.h"
#include "Cubism/Core/Vector.h"
#include "Cubism/Grid/Cartesian.h"
#include <cassert>
#include <mpi.h>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Grid)
/**
 * @defgroup MPI Message Passing Interface
 * @rst
 * The members of this group are parallelized using the distributed memory model
 * implemented by the MPI library.
 * @endrst
 */

/**
 * @ingroup MPI
 * @brief Cartesian MPI block (tensor) field
 * @tparam T Field data type
 * @tparam Mesh Mesh type to be associated with fields
 * @tparam Entity Entity type
 * @tparam RANK Rank of (tensor) fields
 * @tparam Alloc Allocator for field data
 *
 * @rst
 * Cartesian topology composed of block :ref:`field` for the specified entity
 * type.  As opposed to an individual block :ref:`field`, this class manages a
 * structure of arrays (SoA) memory layout for *all* the blocks in the rank
 * local Cartesian topology instead of just individual blocks.  See the
 * :ref:`cartesian` grid section for a non-distributed variant of this class as
 * well as the ``UserState`` extension.
 * @endrst
 */
template <typename T,
          typename Mesh,
          Cubism::EntityType Entity = Cubism::EntityType::Cell,
          size_t RANK = 0,
          typename UserState = Block::FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
class CartesianMPI
    : public Cubism::Grid::Cartesian<T, Mesh, Entity, RANK, UserState, Alloc>
{
    using BaseGrid =
        Cubism::Grid::Cartesian<T, Mesh, Entity, RANK, UserState, Alloc>;
    using IntVec = typename Core::Vector<int, BaseGrid::Dim>;

    using BaseGrid::block_cells_;
    using BaseGrid::block_range_;
    using BaseGrid::global_mesh_;
    using BaseGrid::mesh_;
    using BaseGrid::nblocks_;

public:
    using typename BaseGrid::BaseType;
    using typename BaseGrid::DataType;
    using typename BaseGrid::FieldContainer;
    using typename BaseGrid::FieldState;
    using typename BaseGrid::IndexRangeType;
    using typename BaseGrid::MeshType;
    using typename BaseGrid::MultiIndex;
    using typename BaseGrid::PointType;
    using typename BaseGrid::RangeType;
    using typename BaseGrid::RealType;

    /**
     * @brief Main constructor for a Cartesian MPI block field topology
     * @param comm MPI communicator for this Cartesian topology
     * @param nprocs Number of MPI processes in each dimension
     * @param nblocks Number of blocks per rank
     * @param block_cells Number of cells in each block
     * @param begin Physical origin for the full (all ranks) Cartesian grid
     *              (lower left)
     * @param end Physical end for the full (all ranks) Cartesian grid (top
     *            right)
     * @param gbegin Global begin of physical domain
     * @param gend Global end of physical domain
     */
    CartesianMPI(const MPI_Comm &comm,
                 const MultiIndex &nprocs,
                 const MultiIndex &nblocks,
                 const MultiIndex &block_cells,
                 const PointType &begin = PointType(0),
                 const PointType &end = PointType(1),
                 const PointType &gbegin = PointType(0),
                 const PointType &gend = PointType(1))
        : BaseGrid(), comm_(comm), comm_cart_(MPI_COMM_NULL), nprocs_(nprocs)
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
        const MultiIndex bbegin_rank = rank_index_ * nblocks_;
        block_range_ = IndexRangeType(bbegin_rank, bbegin_rank + nblocks_);

        // mesh and data topology for this rank
        const PointType extent_rank = (end - begin) / PointType(nprocs_);
        const PointType begin_rank =
            begin + PointType(rank_index_) * extent_rank; // rank domain begin
        const PointType end_rank = begin_rank + extent_rank; // rank domain end
        this->initTopology_(gbegin, gend, begin_rank, end_rank, nprocs_);

        // setup global mesh
        const MultiIndex global_blocks = this->getGlobalSize();
        global_mesh_ = new Mesh(mesh_->getGlobalBegin(),
                                mesh_->getGlobalEnd(),
                                block_cells_ * global_blocks,
                                Cubism::MeshIntegrity::FullMesh);
    }

    /** @brief Default constructor */
    CartesianMPI() = default;
    /** @brief Deleted copy constructor */
    CartesianMPI(const CartesianMPI &c) = delete;
    /** @brief Deleted move constructor */
    CartesianMPI(CartesianMPI &&c) = delete;
    /** @brief Default copy assignment */
    CartesianMPI &operator=(const CartesianMPI &c) = default;
    /** @brief Deleted move assignment */
    CartesianMPI &operator=(CartesianMPI &&c) = delete;

    ~CartesianMPI() override
    {
        if (global_mesh_) {
            delete global_mesh_;
            global_mesh_ = nullptr;
        }
    }

    /**
     * @brief Global size of the grid in all dimensions
     * @return Number of blocks in all dimensions in the global grid
     */
    MultiIndex getGlobalSize() const override { return nprocs_ * nblocks_; }
    /**
     * @brief MPI processes in the topology
     * @return Number of processes in the MPI topology
     */
    MultiIndex getNumProcs() const { return nprocs_; }
    /**
     * @brief Cartesian index of MPI process
     * @return Multi-dimensional index of the process
     */
    MultiIndex getProcIndex() const { return rank_index_; }
    /**
     * @brief MPI rank
     * @return Rank of the process in the Cartesian communicator
     */
    int getCartRank() const { return rank_cart_; }
    /**
     * @brief MPI communicator
     * @return Cartesian MPI communicator
     */
    MPI_Comm getCartComm() const { return comm_cart_; }
    /**
     * @brief Test for root process
     * @return True if this is the root rank
     */
    bool isRoot() const { return (0 == rank_cart_); }

private:
    MPI_Comm comm_;         // World communicator
    MPI_Comm comm_cart_;    // Cartesian communicator
    MultiIndex nprocs_;     // Number of MPI processes
    MultiIndex rank_index_; // Cartesian index of this rank
    int rank_cart_;         // Cartesian MPI rank
};

NAMESPACE_END(Grid)
NAMESPACE_END(Cubism)

#endif /* CARTESIANMPI_H_INHL4O2K */
