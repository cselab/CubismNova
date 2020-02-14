// File       : FieldAndLab.h
// Created    : Fri Feb 14 2020 07:34:09 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Convenience data class for testing
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef FIELDANDLAB_H_ZULQKGTJ
#define FIELDANDLAB_H_ZULQKGTJ

#include "Cubism/Block/DataLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"

template <typename T, Cubism::EntityType Entity, size_t DIM>
class FieldAndLab
{
public:
    using Field = Cubism::Block::Field<T, Entity, DIM>;
    using IRange = typename Field::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;
    using DataLab = Cubism::Block::DataLab<Field>;
    using Stencil = typename DataLab::StencilType;
    using DataType = typename Field::DataType;
    using BCVector = typename Field::BCVector;

    enum class Tensorial { Off = 0, On };

    // Changing the stencil extent will break test cases that depend on it.  To
    // simplify tests a symmetric stencil is used.
    FieldAndLab(const Tensorial t = Tensorial::Off)
        : field_(IRange(16)), stencil_(-3, 4, static_cast<bool>(t))
    {
        DataType k = 0;
        for (auto &c : field_) {
            c = k;
            k += 1;
        }

        lab_.allocate(stencil_, field_.getIndexRange());
    }

    Field &getField() { return field_; }
    const Field &getField() const { return field_; }
    DataLab &getLab() { return lab_; }
    const DataLab &getLab() const { return lab_; }
    Stencil &getStencil() { return stencil_; }
    const Stencil &getStencil() const { return stencil_; }

    void loadData(const BCVector *bcs = nullptr)
    {
        auto fields = [this](const MIndex &) -> const Field & {
            return this->field_;
        };
        if (bcs) {
            lab_.loadData(MIndex(0), fields, *bcs);
        } else {
            lab_.loadData(MIndex(0), fields);
        }
    }

private:
    Field field_;
    Stencil stencil_;
    DataLab lab_;
};

#endif /* FIELDANDLAB_H_ZULQKGTJ */
