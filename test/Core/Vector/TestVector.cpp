// File       : TestVector.cpp
// Created    : Mon Apr 22 2019 01:33:29 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Test Vector type
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Core/Timer.h"
#include "Core/Vector.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;
using namespace Cubism;

template <typename T, size_t DIM>
T sumVector(const Vector<T, DIM> &v)
{
    T sum = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += v[i];
    }
    return sum;
}

template <typename T>
T sumArray(const T *p, const size_t n)
{
    T sum = 0;
    for (size_t i = 0; i < n; ++i) {
        sum += *p++;
    }
    return sum;
}

template <typename TA, size_t DIMA, typename TB, size_t DIMB>
void testConstructor()
{
    using VecA = Vector<TA, DIMA>;
    using VecB = Vector<TB, DIMB>;
    using ArrA = typename VecA::ArrayType;
    using ArrB = typename VecB::ArrayType;
    using DataTypeA = typename VecA::DataType;
    using DataTypeB = typename VecB::DataType;

    assert(VecA::Dim > 2);
    assert(VecA::Dim < 4096);

    assert(VecB::Dim > 2);


    cout << "testConstructor():" << '\n';
    cout << '\t' << "typeid(TA).name()       = " << typeid(TA).name() << '\n';
    cout << '\t' << "typeid(TB).name()       = " << typeid(TB).name() << '\n';
    cout << '\t' << "sizeof(VecA)            = " << sizeof(VecA) << '\n';
    cout << '\t' << "sizeof(VecA::ArrayType) = " << sizeof(ArrA) << '\n';
    cout << '\t' << "sizeof(VecB)            = " << sizeof(VecB) << '\n';
    cout << '\t' << "sizeof(VecB::ArrayType) = " << sizeof(ArrB) << '\n';

    vector<VecA> cA; // container A
    vector<VecB> cB; // container B

    Timer t0;
    t0.start();

    // default construction
    VecA a0;
    VecB b0;
    a0[0] = 1;
    b0[0] = 1;

    // copy construction
    cA.push_back(a0);
    cB.push_back(b0);
    cA.back()[0] = 2;
    cB.back()[0] = 2;
    assert(1 == sumVector(a0));
    assert(1 == sumVector(b0));
    assert(2 == sumVector(cA.back()));
    assert(2 == sumVector(cB.back()));

    // move construction
    cA.push_back(std::move(a0));
    cB.push_back(std::move(b0));
    a0[0] = 2; // std::array is copied when moved, still accessible
    b0[0] = 3; // std::array is copied when moved, still accessible
    assert(2 == sumVector(a0));
    assert(3 == sumVector(b0));
    assert(1 == sumVector(cA.back()));
    assert(1 == sumVector(cB.back()));

    // copy construct from underlying array type
    ArrA aa0{0, 1, 2};
    ArrB ab0{1, 2};
    cA.emplace_back(VecA(aa0));
    cB.emplace_back(VecB(ab0));
    assert(3 == sumVector(cA.back()));
    assert(3 == sumVector(cB.back()));

    // move construct from underlying array type
    cA.emplace_back(VecA(std::move(aa0)));
    cB.emplace_back(VecB(std::move(ab0)));
    assert(3 == sumVector(cA.back()));
    assert(3 == sumVector(cB.back()));

    // construct from arbitrary vector type
    cA.push_back(VecA(b0));
    cB.push_back(VecB(a0));
    assert(3 == sumVector(cA.back()));
    assert(2 == sumVector(cB.back()));

    // construct from arbitrary scalar
    cA.push_back(VecA(static_cast<DataTypeA>(1.0)));
    cB.push_back(VecB(static_cast<DataTypeB>(2)));
    assert(VecA::Dim == sumVector(cA.back()));
    assert(2 * VecB::Dim == sumVector(cB.back()));

    // construct from arbitrary initializer list
    cA.push_back(VecA({1.0, 2.0, 3.0}));
    cB.push_back(VecB{3.0f, 4.0f});
    assert(6 == sumVector(cA.back()));
    assert(7 == sumVector(cB.back()));

    // construct from arbitrary std::vector
    vector<int> vi(4096, 1.0);
    vector<float> vf(2, 1.0);
    cA.push_back(VecA(vi));
    cB.push_back(VecB(vf));
    assert(VecA::Dim == sumVector(cA.back()));
    assert(2 == sumVector(cB.back()));

    // construct from arbitrary std::array
    array<size_t, 1> ai{1};
    array<long double, 32> ad{0, 1, 2};
    cA.push_back(VecA(ai));
    cB.push_back(VecB(ad));
    assert(1 == sumVector(cA.back()));
    assert(3 == sumVector(cB.back()));

    // construct from arbitrary pointer type
    vector<char> vc(4096, 1.0);
    vector<unsigned int> vui(1, 1.0);
    cA.push_back(VecA(vc.data(), 4096));
    cB.push_back(VecB(vui.data(), 1));
    assert(VecA::Dim == sumVector(cA.back()));
    assert(1 == sumVector(cB.back()));

    const double telap = t0.stop();

    double res = 0;
    for (size_t i = 0; i < cA.size(); ++i) {
        res += sumVector(cA[i]);
        res += sumVector(cB[i]);
    }

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename TA, size_t DIMA, typename TB, size_t DIMB>
void testAssignment()
{
    using VecA = Vector<TA, DIMA>;
    using VecB = Vector<TB, DIMB>;
    using ArrA = typename VecA::ArrayType;
    using ArrB = typename VecB::ArrayType;

    assert(VecA::Dim > 2);
    assert(VecB::Dim > 2);

    cout << "testAssignment():" << '\n';
    cout << '\t' << "typeid(TA).name()       = " << typeid(TA).name() << '\n';
    cout << '\t' << "typeid(TB).name()       = " << typeid(TB).name() << '\n';
    cout << '\t' << "sizeof(VecA)            = " << sizeof(VecA) << '\n';
    cout << '\t' << "sizeof(VecA::ArrayType) = " << sizeof(ArrA) << '\n';
    cout << '\t' << "sizeof(VecB)            = " << sizeof(VecB) << '\n';
    cout << '\t' << "sizeof(VecB::ArrayType) = " << sizeof(ArrB) << '\n';

    vector<VecA> cA(6); // container A
    vector<VecB> cB(6); // container B

    for (size_t i = 0; i < cA.size(); ++i) {
        assert(0 == sumVector(cA[i]));
    }
    for (size_t i = 0; i < cB.size(); ++i) {
        assert(0 == sumVector(cB[i]));
    }

    Timer t0;
    t0.start();

    // assign vector
    VecA a0{0, 1, 2};
    VecB b0{2, 3};
    cA[0] = a0;
    cB[0] = b0;
    assert(sumVector(cA[0]) == sumVector(a0));
    assert(sumVector(cB[0]) == sumVector(b0));

    // move assign vector
    cA[1] = std::move(a0);
    cB[1] = std::move(b0);
    assert(sumVector(cA[1]) == sumVector(cA[0]));
    assert(sumVector(cB[1]) == sumVector(cB[0]));

    cA[1][0] = 100;
    cB[1][0] = 100;
    assert(sumVector(cA[1]) != sumVector(cA[0]));
    assert(sumVector(cB[1]) != sumVector(cB[0]));

    // assign array
    ArrA aa0{0, 1, 2};
    ArrB ab0{2, 3};
    cA[2] = aa0;
    cB[2] = ab0;
    assert(sumVector(cA[2]) == sumVector(cA[0]));
    assert(sumVector(cB[2]) == sumVector(cB[0]));

    // move assign array
    cA[3] = std::move(aa0);
    cB[3] = std::move(ab0);
    assert(sumVector(cA[3]) == sumVector(cA[0]));
    assert(sumVector(cB[3]) == sumVector(cB[0]));

    cA[3][0] = 100;
    cB[3][0] = 100;
    assert(sumVector(cA[3]) != sumVector(cA[0]));
    assert(sumVector(cB[3]) != sumVector(cB[0]));

    // assign arbitrary vector
    cA[4] = cB[0];
    cB[4] = cA[0];
    assert(sumVector(cA[4]) == sumVector(cB[0]));
    assert(sumVector(cB[4]) == sumVector(cA[0]));

    // assign scalar
    using DataTypeA = typename VecA::DataType;
    using DataTypeB = typename VecB::DataType;
    cA[5] = static_cast<DataTypeA>(2);
    cB[5] = static_cast<DataTypeB>(3);
    assert(sumVector(cA[5]) == 2 * DIMA);
    assert(sumVector(cB[5]) == 3 * DIMB);

    const double telap = t0.stop();

    double res = 0;
    for (size_t i = 0; i < cA.size(); ++i) {
        res += sumVector(cA[i]);
        res += sumVector(cB[i]);
    }

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM>
void testSwap()
{
    using Vec = Vector<T, DIM>;
    using Arr = typename Vec::ArrayType;

    cout << "testSwap():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';
    cout << '\t' << "sizeof(Vec::ArrayType) = " << sizeof(Arr) << '\n';

    vector<Vec> c;

    Timer t0;
    t0.start();

    c.push_back(Vec{1, 0});
    c.push_back(Vec{2, 0});
    assert(1 == sumVector(c[0]));
    assert(2 == sumVector(c[1]));

    c[0].swap(c[1]);
    assert(2 == sumVector(c[0]));
    assert(1 == sumVector(c[1]));

    const double telap = t0.stop();

    double res = 0;
    for (size_t i = 0; i < c.size(); ++i) {
        res += sumVector(c[i]);
    }

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM>
void testLogic()
{
    using Vec = Vector<T, DIM>;
    using Arr = typename Vec::ArrayType;

    assert(Vec::Dim > 1);

    cout << "testLogic():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';
    cout << '\t' << "sizeof(Vec::ArrayType) = " << sizeof(Arr) << '\n';

    Timer t0;
    t0.start();

    Vec v0 = {0};
    Vec v1 = {1, 0};
    Vec v2 = {2, 0};
    Vec v3 = {3, 3};

    // vector logic (Vec::Dim >= 2)
    assert(v0 != v1 && v1 != v2 && v2 != v3);
    assert(v1 != v0 && v2 != v1 && v3 != v2);
    assert(!(v0 == v1) && !(v1 == v2) && !(v2 == v3));
    assert(!(v1 == v0) && !(v2 == v1) && !(v3 == v2));
    assert(v0 <= v1 && v1 <= v2 && v2 <= v3);
    assert(!(v0 >= v1) && !(v1 >= v2) && !(v2 >= v3));
    assert(v3 >= v2 && v2 >= v1 && v1 >= v0);
    assert(!(v3 <= v2) && !(v2 <= v1) && !(v1 <= v0));
    if (2 == Vec::Dim) {
        // Cases when all components are (possibly) different form zero
        assert(v0 < v3 && v1 < v3 && v2 < v3);
        assert(!(v0 > v3) && !(v1 > v3) && !(v2 > v3));
        assert(v3 > v2 && v3 > v1 && v3 > v0);
        assert(!(v3 < v2) && !(v3 < v1) && !(v3 < v0));
    }

    v2[0] = 1;
    assert(v1 == v2);
    assert(v2 == v1);

    // array logic (lexicographic based on Vec::DataType)
    v1[0] = 0;
    v1[1] = 1;
    v2[0] = 0;
    v2[1] = 2;
    const Arr &a0 = v0.getArray();
    const Arr &a1 = v1.getArray();
    const Arr &a2 = v2.getArray();
    const Arr &a3 = v3.getArray();

    assert(a0 != a1 && a1 != a2 && a2 != a3);
    assert(!(a0 == a1) && !(a1 == a2) && !(a2 == a3));
    assert(a0 < a1 && a1 < a2 && a2 < a3);
    assert(a0 <= a1 && a1 <= a2 && a2 <= a3);
    assert(!(a0 > a1) && !(a1 > a2) && !(a2 > a3));
    assert(!(a0 >= a1) && !(a1 >= a2) && !(a2 >= a3));

    const double telap = t0.stop();

    double res = 0;
    res += sumVector(v0);
    res += sumVector(v1);
    res += sumVector(v2);
    res += sumVector(v3);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM>
void testCast()
{
    using Vec = Vector<T, DIM>;
    using Arr = typename Vec::ArrayType;
    using DataType = typename Vec::DataType;

    cout << "testCast():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';
    cout << '\t' << "sizeof(Vec::ArrayType) = " << sizeof(Arr) << '\n';

    Timer t0;
    t0.start();

    Vec v0;

    Arr a0 = static_cast<Arr>(v0);
    const Arr a1 = static_cast<Arr>(v0);
    assert(sumArray(a0.data(), a0.size()) == sumVector(v0));
    assert(sumArray(a1.data(), a1.size()) == sumVector(v0));

    DataType *p0 = static_cast<DataType *>(v0);
    const DataType *p1 = static_cast<const DataType *>(v0);
    assert(sumArray(p0, v0.size()) == sumVector(v0));
    assert(sumArray(p1, v0.size()) == sumVector(v0));

    const double telap = t0.stop();

    double res = 0;
    res += sumVector(v0);
    res += sumArray(a0.data(), a0.size());
    res += sumArray(a1.data(), a1.size());
    res += sumArray(p0, v0.size());
    res += sumArray(p1, v0.size());

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM, bool SUB = true, bool DIV = true>
void testArithmetic()
{
    using Vec = Vector<T, DIM>;
    using DataType = typename Vec::DataType;

    cout << "testArithmetic():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';

    Timer t0;
    t0.start();

    Vec v0;
    Vec v1;
    const size_t n = v0.size();
    const DataType N = static_cast<DataType>(n);

    // both vector operands
    for (size_t i = 0; i < n; ++i) {
        v0[i] = i + 1;
        v1[i] = i + 2;
    }

    // unsigned integral types fail at subtractions from here on
    v0 += v1;
    assert(sumVector(v0) == (N * (N + 1) + N));

    Vec v2 = v0 + v1;
    assert(sumVector(v2) == (3 * (N * (N + 1)) / 2 + 2 * N));

    if (SUB) {
        v2 -= v1;
        assert(sumVector(v2) == (N * (N + 1) + N));

        Vec v3 = v2 - v1;
        assert(sumVector(v3) == ((N * (N + 1)) / 2));
    }

    const DataType fac = 2;

    for (size_t i = 0; i < n; ++i) {
        v0[i] = 1;
    }
    if (SUB) {
        v1 = -v0;
    } else {
        v1 = fac * v0;
    }

    v2 = v0 * v1;
    if (SUB) {
        assert(sumVector(v2) == -N);
    } else {
        assert(sumVector(v2) == fac * N);
    }

    v0 *= v1;
    if (SUB) {
        assert(sumVector(v0) == -N);
    } else {
        assert(sumVector(v2) == fac * N);
    }

    v0 = v1 / v2; // also works for integer division
    assert(sumVector(v0) == N);

    if (SUB) {
        v0 /= v1;
        assert(sumVector(v0) == -N);
    } else {
        v1 /= v0;
        assert(sumVector(v1) == fac * N);
    }

    // rhs scalar operand
    for (size_t i = 0; i < n; ++i) {
        v0[i] = 1;
        v1[i] = 1;
    }

    v0 += fac;
    assert(sumVector(v0) == (fac + 1) * N);

    v0 -= fac; // works for unsigned integral types (only for testing)
    assert(sumVector(v0) == N);

    v0 *= fac;
    assert(sumVector(v0) == fac * N);

    v0 /= fac;
    assert(sumVector(v0) == N);

    v1 = v0 + fac;
    assert(sumVector(v1) == (fac + 1) * N);

    v0 = v1 - fac; // works for unsigned integral types (only for testing)
    assert(sumVector(v0) == N);

    v1 = v0 * fac;
    assert(sumVector(v1) == fac * N);

    v0 = v1 / fac;
    assert(sumVector(v0) == N);

    // lhs scalar operand
    v1 = fac + v0;
    assert(sumVector(v1) == (fac + 1) * N);

    if (SUB) {
        v0 = fac - v1;
        assert(sumVector(v0) == -N);

        v1 = fac * v0;
        assert(sumVector(v1) == -fac * N);

        v0 = fac / v1;
        assert(sumVector(v0) == -N);
    } else {
        v1 = fac - v0;
        assert(sumVector(v1) == N);

        v1 = fac * v0;
        assert(sumVector(v1) == fac * N);

        v0 = fac / v1;
        assert(sumVector(v0) == N);
    }

    if (SUB && DIV) {
        v0 *= -static_cast<DataType>(2);
        v0 = static_cast<DataType>(1) / v0;
        assert(sumVector(v0) == N / static_cast<DataType>(2));
    }

    const double telap = t0.stop();

    double res = 0;
    res += sumVector(v0);
    res += sumVector(v1);
    res += sumVector(v2);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM, bool SUB = true>
void testCommonVecOp()
{
    using Vec = Vector<T, DIM>;
    using DataType = typename Vec::DataType;

    cout << "testCommonVecOp():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';

    Timer t0;
    t0.start();

    // Common vector operation that do not involve real type results

    const DataType nul = static_cast<DataType>(0);
    const DataType one = static_cast<DataType>(1);
    const DataType two = static_cast<DataType>(2);
    const DataType tre = static_cast<DataType>(3);
    Vec v0(two);
    Vec v1(tre);
    const size_t n = v0.size();
    const DataType N = static_cast<DataType>(n);

    assert(v0.normL1() == two * N);
    assert(v0.normLinf() == two);
    assert(v0.normsq() == two * two * N);
    assert(v0.dot(v1) == two * tre * N);
    assert(v0.distsq(v1) == N); // works with unsigned integral for this test
    assert(v0.sum() == two * N);
    assert(v0.prod() == std::pow(two, n));
    assert(v0.min() == two);
    assert(v0.max() == two);
    assert(v0.argmin() == static_cast<size_t>(N - one));
    assert(v0.argmax() == static_cast<size_t>(N - one));

    v0 = Vec{1, 0};
    v1 = Vec{0, 1};
    assert(v0.dot(v1) == nul);

    for (size_t i = 0; i < n; ++i) {
        v0[i] = i + 1;
    }
    assert(v0.min() == one);
    assert(v0.max() == N);
    assert(v0.argmin() == static_cast<size_t>(0));
    assert(v0.argmax() == static_cast<size_t>(N - one));

    assert(sumVector(v0.abs()) == (N * (N + one)) / two);
    if (SUB) {
        v0 = -two;
        assert(sumVector(v0.abs()) == two * N);
    }

    const double telap = t0.stop();

    double res = 0;
    res += sumVector(v0);
    res += sumVector(v1);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM>
void testCommonVecOpReal()
{
    using Vec = Vector<T, DIM>;
    using DataType = typename Vec::DataType;

    cout << "testCommonVecOpReal():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';

    Timer t0;
    t0.start();

    // Common vector operation that do involve real type results

    const DataType one = static_cast<DataType>(1);
    const DataType two = static_cast<DataType>(2);
    const DataType tre = static_cast<DataType>(3);
    Vec v0(two);
    Vec v1(tre);
    Vec vOne(one);

    assert(v0.normL2() == std::sqrt(v0.dot(v0)));
    assert(v0.norm() == std::sqrt(v0.dot(v0))); // alias for normL2
    assert(v0.dist(v1) == vOne.norm());

    const double telap = t0.stop();

    double res = 0;
    res += sumVector(v0);
    res += sumVector(v1);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T>
void testCrossProductVec3()
{
    using Vec = Vector<T, 3>;

    cout << "testCrossProductVec3():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';

    Timer t0;
    t0.start();

    Vec v0{1, 0};
    Vec v1{0, 1};
    Vec v2 = v0.cross(v1);
    assert(v2 == Vec({0, 0, 1}));
    v0[1] = 1;
    v1[0] = 1;
    v2 = v0.cross(v1);
    assert(v2 == Vec({0}));

    const double telap = t0.stop();

    double res = 0;
    res += sumVector(v0);
    res += sumVector(v1);
    res += sumVector(v2);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM>
void testCrossProductVec23()
{
    using Vec = Vector<T, DIM>;
    using DataType = typename Vec::DataType;

    cout << "testCrossProductVec23():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';

    assert(2 == Vec::Dim || 3 == Vec::Dim);

    Timer t0;
    t0.start();

    Vec v0{1, 0};
    Vec v1{0, 1};

    const DataType c3 = v0.getCrossThird(v1);
    assert(c3 == static_cast<DataType>(1));

    const double telap = t0.stop();

    double res = c3;
    res += sumVector(v0);
    res += sumVector(v1);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

template <typename T, size_t DIM>
void testIterators()
{
    using Vec = Vector<T, DIM>;
    using DataType = typename Vec::DataType;

    cout << "testIterators():" << '\n';
    cout << '\t' << "typeid(T).name()       = " << typeid(T).name() << '\n';
    cout << '\t' << "sizeof(Vec)            = " << sizeof(Vec) << '\n';

    Timer t0;
    t0.start();

    Vec v0;

    // forward
    for (DataType &v : v0) {
        v = 1;
    }
    assert(sumVector(v0) == v0.size());

    DataType sum = 0;
    for (const DataType &v : v0) {
        // v += 1; // compile-time error
        sum += v;
    }
    assert(sumVector(v0) == sum);

    // reverse
    for (auto v = v0.rbegin(); v < v0.rend(); ++v) {
        *v += 1;
    }
    assert(sumVector(v0) == 2 * v0.size());

    for (auto v = v0.crbegin(); v < v0.crend(); ++v) {
        // *v += 1; // compile-time error
        sum += *v;
    }
    assert(sumVector(v0) == 2 * v0.size());
    assert(sum == 3 * v0.size());

    const double telap = t0.stop();

    double res = sum;
    res += sumVector(v0);

    cout << "Result: " << res << '\n';
    cout << "Took:   " << telap << " seconds" << '\n';
}

int main(void)
{
    // construction
    testConstructor<double, 3, int, 10>();
    testConstructor<ptrdiff_t, 10, float, 3>();

    // assignment
    testAssignment<double, 3, int, 10>();
    testAssignment<ptrdiff_t, 10, float, 3>();

    // swap
    testSwap<double, 4>();
    testSwap<size_t, 16>();

    // logic operators
    testLogic<double, 2>();
    testLogic<float, 3>();
    testLogic<ptrdiff_t, 5>();

    // casts
    testCast<long double, 3>();
    testCast<int, 1>();

    // aritmetic
    testArithmetic<float, 3>();
    testArithmetic<double, 10>();
    testArithmetic<long double, 10>();
    testArithmetic<int, 16, true, false>();    // no fractional integer
                                               // division
                                               // test
    testArithmetic<size_t, 8, false, false>(); // no fractional integer
                                               // division test and no
                                               // subtraction tests

    // common vector operations
    testCommonVecOp<float, 2>();
    testCommonVecOp<double, 16>();
    testCommonVecOp<long double, 9>();
    testCommonVecOp<int, 6>();
    testCommonVecOp<unsigned int, 6, false>();

    // common vector operations with fractional return type
    testCommonVecOpReal<float, 23>();
    testCommonVecOpReal<double, 101>();
    testCommonVecOpReal<long double, 3>();
    // testCommonVecOpReal<long, 3>(); // link-time error

    // test cross product
    testCrossProductVec23<float, 2>();
    testCrossProductVec23<double, 3>();
    testCrossProductVec23<int, 3>();
    // testCrossProductVec23<double, 4>(); // compile-time error

    testCrossProductVec3<float>();
    testCrossProductVec3<double>();
    testCrossProductVec3<int>();
    testCrossProductVec3<unsigned int>(); // ok for positive number results

    // test iterators
    testIterators<float, 4>();
    testIterators<size_t, 9>();

    return 0;
}
