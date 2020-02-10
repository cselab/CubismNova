// File       : DataLab.h
// Created    : Mon Feb 10 2020 06:53:39 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Data laboratory with stencil specification
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef DATALAB_H_O091Y6A2
#define DATALAB_H_O091Y6A2

#include "Cubism/Block/Data.h"

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename BlockAlloc = AlignedBlockAllocator<T>>
class DataLab : public Data<T, Entity, DIM, BlockAlloc>
{
public:
    DataLab() {}
    ~DataLab() override {}

    // stencil struct
    void setStencil(/* Stencil */);
    // allocates the memory
    // can not call this method if data loaded flag is set

    void loadData(/* Field */);
    // fills the memory

    void clear();
    // flag that data can be destroyed

private:
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* DATALAB_H_O091Y6A2 */
