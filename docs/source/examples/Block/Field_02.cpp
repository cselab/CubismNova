#include <Cubism/Block/Field.h>

int main(void)
{
    // CUBISM_DIMENSION-ional field
    using CellField = Cubism::Block::CellField<float>;
    using FieldView = Cubism::Block::FieldView<CellField>; // a CellField view
    using IRange = typename CellField::IndexRangeType;
    using MIndex = typename IRange::MultiIndex;

    MIndex cells(16);          // number of cells in field (256 cells)
    IRange cell_domain(cells); // index range spanned by cell domain
    CellField cf(cell_domain); // cell field (memory untouched, owner)
    CellField co(cell_domain); // some other cell field, memory owner

    FieldView cv(cf); // cv views into cf
    cv.setView(co);   // cv now views into co (cheap)
    cv.copyData(cf);  // data of cf is copied into co (expensive)

    CellField cc = cv.copy(); // create full copy of co in cc (move assignment)
    FieldView ccv(cc);        // another view
    cf = ccv;                 // assign to field from view

    return 0;
}
