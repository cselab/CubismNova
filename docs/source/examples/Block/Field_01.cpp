#include "Block/Field.h"

using namespace Cubism;

int main(void)
{
    // a custom field state data structure
    struct MyFieldState {
        size_t rank; // required
        size_t comp; // required
        size_t myID; // my addition
    };

    using CellField = Block::CellField<double, MyFieldState, 2>; // 2D field
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
