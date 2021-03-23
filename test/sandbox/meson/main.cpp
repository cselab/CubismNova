// File       : main.cpp
// Created    : Tue Mar 23 2021 09:51:59 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Test program using the CubismNova library
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include <Cubism/Block/Field.h>

int main(void)
{
    // a custom field state data structure
    struct MyFieldState {
        size_t myID; // custom field state
    };

    constexpr size_t dim = 2; // 2D problem
    using CellField = Cubism::Block::CellField<double, dim, MyFieldState>;
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);          // number of cells in field (256 cells)
    IRange cell_domain(cells); // index range spanned by cell domain
    CellField cf(cell_domain); // cell field (memory untouched)

    cf.getState().myID = 101; // set my custom state
    for (auto &c : cf) {      // cell c in cell field cf
        c = 0;
    }

    return 0;
}
