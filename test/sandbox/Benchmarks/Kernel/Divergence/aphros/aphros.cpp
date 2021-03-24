// File       : aphros.cpp
// Created    : Wed Mar 04 2020 04:35:07 PM (-0800)
// Author     : Fabian Wermelinger
// Description: Compute divergence test kernel
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include <cmath>
#include <cstdio>
#include <sstream>

#include "distr/distrsolver.h"
#include "kernel/hydro.h"
#include "kernel/kernelmeshpar.h"

#include "../../../Utils/Timer.h"

using M = MeshStructured<double, 3>;
using K = Hydro<M>;
using Par = typename K::Par;

template <typename M_, typename K_>
class MySolver : public DistrSolver<M_, K_>
{
public:
    MySolver(MPI_Comm comm, Vars &var0, Par &par) {}
    ~MySolver() {}

private:
};

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    std::string myconf = R"EOF(
set int bx 8
set int by 8
set int bz 8

set int bsx 32
set int bsy 32
set int bsz 32

set int px 1
set int py 1
set int pz 1

set int CHECKNAN 0
set int dim 3

set string dumpformat default

set int comm_size 1

set string backend cubismnc
set int histogram 0
set int openmp 0
set int mpi_compress_msg 0

set int max_step 1
set int num_frames 1

set int hl 2
set int hypre_print 0
set double hypre_symm_tol 1e-12
set double hypre_vort_tol 1e-12
set double hypre_gen_tol 1e-12
set int periodic 1
set int hypre_periodic_x 1
set int hypre_periodic_y 1
set int hypre_periodic_z 1
set string hypre_gen_solver gmres

set int loc_periodic_x 1
set int loc_periodic_y 1
set int loc_periodic_z 1

set int verbose 0
set int output 0
set int verbose_stages 0
set int verbose_time 0
set int verbose_openmp 0

set int iter 1

set double extent 1.

set string hypre_symm_solver pcg
set string hypre_vort_solver pcg
set int hypre_symm_maxiter 100
set int hypre_vort_maxiter 100
set int hypre_gen_maxiter 30
)EOF";
    std::stringstream conf;
    conf << "\n" << myconf;

    Vars var;       // parameter storage
    Parser ip(var); // parser
    ip.RunAll(conf);
    Par par;
    DistrSolver<M, K> ds(MPI_COMM_WORLD, var, par);

    MPI_Finalize();
    return 0;
}
