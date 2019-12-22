// Examples of client code for possible support

class M;             // mesh
using Scal = double; // scalar
class Vect;          // vector
class IdxCell;       // cell index
class IdxFace;       // face index
template <class T, class Idx>   // field
class Field;


Scal getTotalVolume(const M& m) {
    Scal s(0);
    for (auto c : m.Cells()) {
        s += m.getVolume(c);
    }
    return s;
}

Vect getFirstMoment(const Field<Scal, IdxCell>& fc, const M& m) {
    Vect x(0);
    for (auto c : m.Cells()) {
        x += m.getCenter(x) * (fc[c] * m.getVolume(x));
    }
    return x;
}

void interpolate(const Field<Scal, IdxFace>& ff,
                 Field<Scal, IdxCell>& fc,
                 const M& m) {
    fc.Init(m);
    for (auto c : m.Cells()) {
        fc[c] = 0;
        for (size_t i = 0; i < m.getNumFaces(c); ++i) {
            IdxFace f = m.getFace(c, i);
            fc[c] += ff[f];
        }
        fc[c] /= m.getNumFaces(c);
    }
}

void interpolate(const Field<Scal, IdxCell>& fc,
                 Field<Scal, IdxFace>& ff,
                 const M& m) {
    ff.Init(m);
    for (auto f : m.Faces()) {
        ff[f] = 0;
        for (size_t i = 0; i < m.getNumCells(f); ++i) {
            IdxCell c = m.getCell(f, i);
            ff[f] += fc[c];
        }
        ff[f] /= m.getNumCells(f);
    }
}
