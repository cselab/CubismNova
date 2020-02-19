// File       : BlockProcessor.h
// Date       : Fri 01 Apr 2016 05:52:01 PM CEST
// Author     : Fabian Wermelinger
// Description: Process all blocks
// Copyright 2016 ETH Zurich. All Rights Reserved.
#ifndef BLOCKPROCESSORMPI_H_IKFSZWUJ
#define BLOCKPROCESSORMPI_H_IKFSZWUJ

#include <vector>
#ifdef __OPENMP
#include <omp.h>
#endif /* __OPENMP */

#include "Types.h"

#include <iostream>
using namespace std;

template <typename TLab, typename Operator, typename TGrid>
inline void BlockProcessor(Operator rhs,
                           TGrid &grid,
                           const Real t = 0,
                           const bool record = false)
{
    vector<BlockInfo> avail0;

#ifdef __OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif /* __OPENMP */

    TLab * labs = new TLab[nthreads];

    // Setup the static stencil information for this kernel (operator)
    const int ss[3] = {rhs.stencil.sx, rhs.stencil.sy, rhs.stencil.sz};
    const int se[3] = {rhs.stencil.ex, rhs.stencil.ey, rhs.stencil.ez};
    for (int i = 0; i < nthreads; ++i)
        labs[i].prepare(grid, ss, se, rhs.stencil.tensorial);

    // process inner blocks
    avail0 = grid.getBlocksInfo();
    BlockInfo * ary0 = &avail0.front();

#pragma omp parallel num_threads(nthreads)
    {
#ifdef __OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif /* __OPENMP */
        TLab& mylab = labs[tid];

#pragma omp for schedule(dynamic,1)
        for (size_t i = 0; i < avail0.size(); i++) {
            mylab.load(ary0[i], t);
            rhs(mylab, ary0[i], *(FluidBlock*)ary0[i].ptrBlock);
        }
    }

    // clean up
    if(labs!=NULL)
    {
        delete[] labs;
        labs=NULL;
    }
}

#endif /* BLOCKPROCESSORMPI_H_IKFSZWUJ */
