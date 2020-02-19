// File       : test-cubismo.cpp
// Date       : Fri 01 Apr 2016 05:50:11 PM CEST
// Author     : Fabian Wermelinger
// Description: Short Cubism MPI test
// Copyright 2016 ETH Zurich. All Rights Reserved.

#include <cmath>

#include "BlockProcessor.h"
#include "Cubism/StencilInfo.h"
#include "Types.h"
#ifdef _DUMP_
#include "Cubism/HDF5Dumper.h"
#endif /* _DUMP_ */
#include "../../../Utils/Timer.h"
#include <cstdio>

using namespace std;
using namespace cubism;
using Utils::Timer;

const double twopi = atan(1.0)*8.0;

void initGrid(NodeGrid &grid)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();

#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    double pos[3];
                    info.pos(pos, ix, iy, iz);
                    const double x = pos[0];
                    const double y = pos[1];
                    const double z = pos[2];
                    const double q1 = sin(twopi*x) * cos(twopi*y) * sin(twopi*z);
                    const double q2 = cos(twopi*y) * cos(twopi*x) * sin(twopi*z);
                    const double q3 = cos(twopi*z) * cos(twopi*y) * sin(twopi*x);
                    b(ix, iy, iz).data[0] = q1;
                    b(ix, iy, iz).data[1] = q2;
                    b(ix, iy, iz).data[2] = q3;
                    // b(ix, iy, iz).data[3] = 5.0*q1;
                    // b(ix, iy, iz).data[4] = 5.0*q2;
                    // b(ix, iy, iz).data[5] = 5.0*q3;
                    // b(ix, iy, iz).data[6] = 1.0;
                    // b(ix, iy, iz).data[7] = 0.0;
                }
    }
}

struct Evaluate123Divergence_CPP
{
    // second order FD divergence of data[0], data[1], data[2] (of
    // FluidElement)

    StencilInfo stencil;

    // stencil: first 3 are stencil to the left from zero (inclusive)
    //          next 3 are stencil to the right from zero (exclusive)
    //          next 1 says tensorial true/false
    //          next 1 specifies the n fields you need to evaluate the operator
    //          next n are indices of the fields to be communicated by MPI
    Evaluate123Divergence_CPP(): stencil(-1,-1,-1,2,2,2, false, 3, 0,1,2) {}

    // operator interface
    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        typedef BlockType B;
        const Real ih = 1. / info.h_gridpoint;
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const Real dd0_dx = lab(ix+1,iy,iz).data[0] - lab(ix-1,iy,iz).data[0];
                    const Real dd1_dy = lab(ix,iy+1,iz).data[1] - lab(ix,iy-1,iz).data[1];
                    const Real dd2_dz = lab(ix,iy,iz+1).data[2] - lab(ix,iy,iz-1).data[2];

                    // we stuff it into data[3]
                    o(ix, iy, iz).data[3] =
                        0.5 * (dd0_dx * ih + dd1_dy * ih + dd2_dz * ih);
                }
    }
};


template <int _ID>
struct StreamerPickOne_HDF5
{
    static const int NCHANNELS = 1; // we dump scalar fields with this streamer

    // write
    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        output[0] = el.data[_ID];
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        return el.data[_ID];
    }

    static const char * getAttributeName() { return "Scalar"; }
};


int main(int argc, char* argv[])
{
    int provided;
    const int bpdx = 16;
    const int bpdy = bpdx;
    const int bpdz = bpdx;
    const double maxextent = 1.0;

    // allocate the grid
    Timer t;
    NodeGrid *const mygrid = new NodeGrid(bpdx, bpdy, bpdz, maxextent);

    // initialize values
    initGrid(*mygrid);
    const double t0 = t.stop();
#ifdef _DUMP_
    t.start();
    DumpHDF5<StreamerPickOne_HDF5<0>, Real>(*mygrid, 0, 0, "data0");
    DumpHDF5<StreamerPickOne_HDF5<1>, Real>(*mygrid, 0, 0, "data1");
    DumpHDF5<StreamerPickOne_HDF5<2>, Real>(*mygrid, 0, 0, "data2");
    double td = t.stop();
#endif /* _DUMP_ */

    // evaluate div([data[0], data[1], data[2]]')
    Evaluate123Divergence_CPP myOperator;
    t.start();
    BlockProcessor<PLab>(myOperator, *mygrid); // absorbing BC
    const double t1 = t.stop();

    // dump result
#ifdef _DUMP_
    t.start();
    DumpHDF5<StreamerPickOne_HDF5<3>, Real>(*mygrid, 0, 0, "div123");
    td += t.stop();
    printf("%e\t%e\t%e\t%e\n", t0, t1, td, t0 + t1 + td);
#else
    printf("%e\t%e\t%e\n", t0, t1, t0 + t1);
#endif /* _DUMP_ */
    return 0;
}
