// File       : Types.h
// Date       : Fri 01 Apr 2016 05:54:03 PM CEST
// Author     : Fabian Wermelinger
// Description: Some required types
// Copyright 2016 ETH Zurich. All Rights Reserved.
#ifndef TYPES_H_NTYWLJPY
#define TYPES_H_NTYWLJPY

#include "Cubism/BlockInfo.h"
#include "Cubism/BlockLab.h"
#include "Cubism/Grid.h"

#include "BoundaryConditions.h"

using namespace std;
using namespace cubism;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

struct FluidElement
{
    typedef Real RealType;
    Real data[4]; // just an example, could also specify multiple members with
                  // more descriptive names
    FluidElement() {}
    void clear()
    {
        for (int i = 0; i < 4; ++i)
            data[i] = 0.0;
    }
};

struct FluidBlock
{
    static const int sizeX = _BLOCKSIZEX_;
    static const int sizeY = _BLOCKSIZEY_;
    static const int sizeZ = _BLOCKSIZEZ_;

    static const int gptfloats = sizeof(FluidElement)/sizeof(Real);

    typedef typename FluidElement::RealType RealType;
    typedef FluidElement ElementType;
    typedef FluidElement element_type;

    FluidElement __attribute__((__aligned__(CUBISM_ALIGNMENT))) data[_BLOCKSIZEZ_][_BLOCKSIZEY_][_BLOCKSIZEX_];

    Real __attribute__((__aligned__(CUBISM_ALIGNMENT))) tmp[_BLOCKSIZEZ_][_BLOCKSIZEY_][_BLOCKSIZEX_][gptfloats];

    void clear_data()
    {
        const int N = sizeX*sizeY*sizeZ;
        FluidElement * const e = &data[0][0][0];
        for(int i=0; i<N; ++i) e[i].clear();
    }

    void clear_tmp()
    {
        const int N = sizeX * sizeY * sizeZ * gptfloats;

        Real * const e = &tmp[0][0][0][0];
        for(int i=0; i<N; ++i) e[i] = 0;
    }

    void clear()
    {
        clear_data();
        clear_tmp();
    }

    inline FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }

    inline const FluidElement& operator()(int ix, int iy=0, int iz=0) const
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class MyLabAbsorbing: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    MyLabAbsorbing(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "MyTestLab"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        // "this->" references attributes inherited from BlockLab
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};

typedef Grid<FluidBlock, std::allocator> NodeGrid;
typedef MyLabAbsorbing<FluidBlock> ALab; // absorbing BC (Lab)
typedef BlockLab<FluidBlock> PLab;       // periodic BC (Lab)

#endif /* TYPES_H_NTYWLJPY */
