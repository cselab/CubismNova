// File       : TestVector.cpp
// Created    : Mon Apr 22 2019 01:33:29 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Test Vector type
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Core/Timer.h"
#include "Core/Vector.h"

#include <array>
#include <cassert>
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

template <typename TA, size_t DIMA, typename TB, size_t DIMB>
double testConstructor()
{
    using VecA = Vector<TA, DIMA>;
    using VecB = Vector<TB, DIMB>;

    double res = 0.;
    { // default construction
        VecA v0;
        const TA s0 = sumVector(v0);
        assert(s0 == 0);
        res += s0;

        VecB v1;
        const TB s1 = sumVector(v1);
        assert(s1 == 0);
        res += s1;
    }

    { // copy construction
        VecA v0;
        v0[0] = 1;
        VecA v1(v0);
        const TA s1 = sumVector(v1);
        assert(s1 == 1);
        res += s1;
    }

    { // move construction
        VecA v0;
        v0[0] = 1;
        VecA v1(std::move(v0));
        const TA s1 = sumVector(v1);
        assert(v0.data() == nullptr);
        // v0[0] = 1; // segmentation fault
        assert(s1 == 1);
        res += s1;
    }

    { // construct from hash
        VecA v0;
        v0[0] = 1;
        const typename VecA::HashType h = v0.getHash();
        VecA v1(h);
        const TA s1 = sumVector(v1);
        assert(s1 == 1);
        res += s1;
    }

    { // construct from different vector
        VecA v0;
        v0[0] = 1;
        VecB v1(v0);
        const TB s1 = sumVector(v1);
        assert(s1 == 1);
        res += s1;
    }

    { // scalar construction
        VecA v0(static_cast<TA>(1));
        const TA s0 = sumVector(v0);
        assert(s0 == 1 * VecA::Dim);
        res += s0;

        VecB v1(static_cast<TB>(2));
        const TB s1 = sumVector(v1);
        assert(s1 == 2 * VecB::Dim);
        res += s1;
    }

    { // list initializer
        VecA v0{0, 1, 2};
        const TA s0 = sumVector(v0);
        assert(s0 == 3);
        res += s0;

        VecB v1({0, 1, 2});
        const TB s1 = sumVector(v1);
        assert(s1 == 3);
        res += s1;

        VecA v2 = {1, 2, 3};
        const TA s2 = sumVector(v2);
        assert(s2 == 6);
        res += s2;
    }

    { // construct from std::vector
        vector<TB> v0(100, 1);
        VecA v1(v0);
        const TA s1 = sumVector(v1);
        assert(s1 == 1 * VecA::Dim);
        res += s1;
    }

    { // construct from std::array
        array<TB, 2> v0 = {100, 1};
        VecA v1(v0);
        const TA s1 = sumVector(v1);
        assert(s1 == 101);
        res += s1;
    }

    return res;
}

template <typename TA, size_t DIMA, typename TB, size_t DIMB>
double testAssignment()
{
    using VecA = Vector<TA, DIMA>;
    using VecB = Vector<TB, DIMB>;

    double res = 0.;
    { // default assignment
        VecA v0;
        v0[0] = 1;
        const TA s0 = sumVector(v0);
        assert(s0 == 1);
        res += s0;

        VecA v1;
        v1[0] = 7;
        const TA s1 = sumVector(v1);
        assert(s1 == 7);
        res += s1;

        v0 = v1;
        const TA s2 = sumVector(v0);
        assert(s2 == 7);
        res += s2;
    }

    { // default move assignment
        VecA v0;
        v0[0] = 1;
        const TA s0 = sumVector(v0);
        assert(s0 == 1);
        res += s0;

        VecA v1;
        v1[0] = 7;
        const TA s1 = sumVector(v1);
        assert(s1 == 7);
        res += s1;

        v0 = std::move(v1);
        const TA s2 = sumVector(v0);
        assert(v1.data() == nullptr);
        // v1[0] = 1; // segmentation fault
        assert(s2 == 7);
        res += s2;
    }

    { // assignment for different vector types
        VecA v0;
        v0[0] = 1;
        const TA s0 = sumVector(v0);
        assert(s0 == 1);
        res += s0;

        VecB v1;
        v1[0] = 7;
        const TB s1 = sumVector(v1);
        assert(s1 == 7);
        res += s1;

        v0 = v1;
        const TA s2 = sumVector(v0);
        assert(s2 == 7);
        res += s2;
    }

    return res;
}

template <typename T, size_t DIM>
double testSwap()
{
    // swapping only works with identical vector types
    using Vec = Vector<T, DIM>;

    double res = 0.;
    {
        Vec v0;
        // const void *p0 = reinterpret_cast<const void *>(v0.getHash().data());
        v0[0] = 1;
        const T s0 = sumVector(v0);
        assert(s0 == 1);
        res += s0;

        // vector<double> vec0(10, 0);
        // vector<double> vec1(10, 1);
        // std::cout << (void *)vec0.data() << std::endl;
        // std::cout << (void *)vec1.data() << std::endl;
        // swap(vec0, vec1);
        // std::cout << (void *)vec0.data() << std::endl;
        // std::cout << (void *)vec1.data() << std::endl;

        Vec v1;
        const T s1 = sumVector(v1);
        assert(s1 == 0);
        res += s1;

        v1.swap(v0);
        // const void *p1 = reinterpret_cast<const void *>(v1.getHash().data());
        // assert(p1 == p0);
        const T s2 = sumVector(v1);
        assert(s2 == 1);
        res += s2;

        swap(v0, v1);
        const T s3 = sumVector(v1);
        assert(s3 == 0);
        res += s3;
    }

    return res;
}

double testHash()
{
    constexpr size_t DIM = 5;
    using T = int;
    using Vec = Vector<T, DIM>;

    Vec v0({~0, 1, 15});

    typename Vec::HashType h = v0.getHash();
    const string hashStr = Vec::getHashString(h);

    // this test is meaningful only with an integral type
    ostringstream os;
    os << "0x";
    for (size_t i = 0; i < sizeof(T) * (static_cast<int>(DIM) - 3); ++i) {
        os << std::hex << std::setfill('0') << std::setw(2) << 0;
    }
    for (int i = 0; i < static_cast<int>(sizeof(T)) - 1; ++i) {
        os << std::hex << std::setfill('0') << std::setw(2) << 0;
    }
    os << "0f";

    for (int i = 0; i < static_cast<int>(sizeof(T)) - 1; ++i) {
        os << std::hex << std::setfill('0') << std::setw(2) << 0;
    }
    os << "01";

    for (size_t i = 0; i < sizeof(T); ++i) {
        os << "ff";
    }
    assert(os.str() == hashStr);

    cout << '\t' << "HashString = " << hashStr << '\n';

    for (size_t i = 0; i < v0.size(); ++i) {
        cout << '\t' << "v0[" << i << "] = " << Vec::getComponent(i, h) << '\n';
    }

    return sumVector(v0);
}

int main(void)
{
    Timer t0;

    // construction
    {
        t0.start();
        const double res = testConstructor<double, 3, int, 10>();
        const double telapsed = t0.stop();
        cout << std::setw(60) << std::left
             << "testConstructor<double, 3, int, 10>():";
        cout << "result = " << res << "; took " << telapsed << " sec" << '\n';
    }
    {
        t0.start();
        const double res = testConstructor<ptrdiff_t, 10, float, 3>();
        const double telapsed = t0.stop();
        cout << std::setw(60) << std::left
             << "testConstructor<ptrdiff_t, 10, float, 3>():";
        cout << "result = " << res << "; took " << telapsed << " sec" << '\n';
    }

    // assignment
    {
        t0.start();
        const double res = testAssignment<double, 3, int, 10>();
        const double telapsed = t0.stop();
        cout << std::setw(60) << std::left
             << "testAssignment<double, 3, int, 10>():";
        cout << "result = " << res << "; took " << telapsed << " sec" << '\n';
    }
    {
        t0.start();
        const double res = testAssignment<ptrdiff_t, 10, float, 3>();
        const double telapsed = t0.stop();
        cout << std::setw(60) << std::left
             << "testAssignment<ptrdiff_t, 10, float, 3>():";
        cout << "result = " << res << "; took " << telapsed << " sec" << '\n';
    }

    // swap
    {
        t0.start();
        const double res = testSwap<double, 3>();
        const double telapsed = t0.stop();
        cout << std::setw(60) << std::left << "testSwap<double, 3>():";
        cout << "result = " << res << "; took " << telapsed << " sec" << '\n';
    }

    // hash
    {
        t0.start();
        const double res = testHash();
        const double telapsed = t0.stop();
        cout << std::setw(60) << std::left << "testHash():";
        cout << "result = " << res << "; took " << telapsed << " sec" << '\n';
    }

    return 0;
}
