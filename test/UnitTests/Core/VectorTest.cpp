// File       : VectorTest.cpp
// Created    : Tue Apr 30 2019 05:06:37 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Unit tests for Cubism/Core/Vector.h
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Core/Vector.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <vector>

namespace
{
// some test types
template <size_t D>
using V_f = Cubism::Core::Vector<float, D>;
template <size_t D>
using V_d = Cubism::Core::Vector<double, D>;
template <size_t D>
using V_i = Cubism::Core::Vector<int, D>;
template <size_t D>
using V_m = Cubism::Core::Vector<size_t, D>;

template <typename T, size_t DIM>
T sumVector(const Cubism::Core::Vector<T, DIM> &v)
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

// Construction of vectors
TEST(Vector, Construction)
{
    // fundamental constructors
    {
        using Vec3 = V_f<3>;
        using Arr3 = typename Vec3::ArrayType;

        EXPECT_EQ(12, Vec3::Bytes);
        EXPECT_EQ(3, Vec3::Dim);

        Vec3 v0; // default construction
        EXPECT_EQ(3, v0.size());
        EXPECT_EQ(0, sumVector(v0));

        Vec3 v1(v0); // copy construction
        EXPECT_EQ(0, sumVector(v1));

        Vec3 v2(std::move(v0)); // move construction
        v0[0] = 1;
        EXPECT_EQ(1, sumVector(v0)); // std::array does not move the pointers.
                                     // v0 is still accessible and uses
                                     // different memory than v2
        EXPECT_EQ(0, sumVector(v2));

        Arr3 a3{0, 1, 2};
        Vec3 v3(a3); // copy construct from std::array type
        EXPECT_EQ(3, sumVector(v3));

        Vec3 v4(Arr3{1, 2, 3}); // move construct from std::array type
        EXPECT_EQ(6, sumVector(v4));
    }

    // convenience constructors
    {
        // initializer lists
        V_d<9> v0{0, 1, 2, 3}; // remaining elements are initialized to zero
        V_i<2> v1{1, 2};
        V_m<2> v2{1, 2, 3}; // only first two elements are taken into account
        EXPECT_EQ(9, v0.size());
        EXPECT_EQ(2, v1.size());
        EXPECT_EQ(2, v2.size());
        EXPECT_EQ(6, sumVector(v0));
        EXPECT_EQ(3, sumVector(v1));
        EXPECT_EQ(3, sumVector(v2));

        // from arbitrary scalar types
        V_i<6> v3(static_cast<double>(1));
        V_f<3> v4(static_cast<size_t>(2));
        EXPECT_EQ(6, sumVector(v3));
        EXPECT_EQ(6, sumVector(v4));

        // from arbitrary Cubism::Core::Vector types
        V_d<4> v5{2, 3};
        V_i<9> v6(v5); // v6.size() > v5.size()
        V_m<1> v7(v5); // v7.szie() < v5.size()
        EXPECT_EQ(9, v6.size());
        EXPECT_EQ(1, v7.size());
        EXPECT_EQ(5, sumVector(v6));
        EXPECT_EQ(2, sumVector(v7));

        // from arbitrary std::vector types
        std::vector<int> vec0(16, 1);
        V_d<5> v8(vec0);
        V_m<32> v9(vec0);
        EXPECT_EQ(5, sumVector(v8));
        EXPECT_EQ(16, sumVector(v9));

        // from arbitrary std::array types
        std::array<char, 2> arr0{1, 1};
        V_f<1> v10(arr0);
        V_i<3> v11(arr0);
        EXPECT_EQ(1, sumVector(v10));
        EXPECT_EQ(2, sumVector(v11));

        // from arbitrary pointer types
        V_m<1> v12(vec0.data(), vec0.size());
        V_i<3> v13(arr0.data(), arr0.size());
        EXPECT_EQ(1, sumVector(v12));
        EXPECT_EQ(2, sumVector(v13));
    }
}

// Assignment of vectors
TEST(Vector, Assignment)
{
    using Vec3 = V_i<3>;
    using Arr3 = typename Vec3::ArrayType;

    Vec3 zero;
    Vec3 v0;
    Vec3 v1{1, 1, 1};
    EXPECT_EQ(0, sumVector(v0));
    EXPECT_EQ(3, sumVector(v1));

    v0 = v1; // copy assignment
    EXPECT_EQ(3, sumVector(v0));

    Vec3 v2 = Vec3{1, 1, 1}; // move assignment
    EXPECT_EQ(3, sumVector(v2));

    v0 = zero;
    Arr3 a3{1, 1, 1};
    v0 = a3; // copy assignment form std::array
    EXPECT_EQ(3, sumVector(v0));

    v0 = zero;
    v0 = Arr3{1, 1, 1}; // move copy assignment
    EXPECT_EQ(3, sumVector(v0));

    // scalar assignment
    v0 = zero;
    v0 = static_cast<Vec3::DataType>(1);
    EXPECT_EQ(3, sumVector(v0));

    v0 = zero;
    v0 = static_cast<double>(1);
    EXPECT_EQ(3, sumVector(v0));

    // arbitrary vector assignment
    V_f<32> vf32{1};
    v0 = vf32; // if rhs.size() >= v0.size() -- all elements from rhs are
               // assigned to v0
    EXPECT_EQ(1, sumVector(v0));

    V_d<1> vd1{5};
    v1 = vd1; // if rhs.size() < v1.size() -- only rhs.size() elements are
              // copied to v1, the remaining elements are left untouched.
    EXPECT_EQ(7, sumVector(v1)); // 5 + 1 + 1
}

// Swap vectors
TEST(Vector, Swap)
{
    V_i<3> v0;
    V_i<3> v1{1, 1, 1};
    EXPECT_EQ(0, sumVector(v0));
    EXPECT_EQ(3, sumVector(v1));

    v0.swap(v1);
    EXPECT_EQ(3, sumVector(v0));
    EXPECT_EQ(0, sumVector(v1));
}

// Logic
TEST(Vector, Logic)
{
    using Vec3 = V_i<3>; // test with integers
    using Vec2 = V_i<2>;

    // Cubism::Core::Vector logic operators do not follow lexicographical order
    Vec3 v3_0 = {0};
    Vec3 v3_1 = {1, 0};
    Vec3 v3_2 = {2, 0};
    Vec3 v3_3 = {3, 3};
    Vec3 v3_4 = {3, 3};

    EXPECT_NE(v3_0, v3_1);
    EXPECT_NE(v3_1, v3_2);
    EXPECT_NE(v3_2, v3_3);

    // less-than and greater-than operators are expected to fail if there are
    // vector elements that are equal
    EXPECT_FALSE(v3_1 < v3_2);
    EXPECT_FALSE(v3_2 > v3_1);

    EXPECT_LE(v3_0, v3_1);
    EXPECT_LE(v3_1, v3_2);
    EXPECT_LE(v3_2, v3_3);
    EXPECT_LE(v3_3, v3_4);

    EXPECT_GE(v3_1, v3_0);
    EXPECT_GE(v3_2, v3_1);
    EXPECT_GE(v3_3, v3_2);
    EXPECT_GE(v3_3, v3_4);

    Vec2 v2_0 = v3_0;
    Vec2 v2_1 = v3_1;
    Vec2 v2_2 = v3_2;
    Vec2 v2_3 = v3_3;
    Vec2 v2_4 = v3_4;
    EXPECT_LT(v2_0, v2_3); // less-than (all elements must be less-than)
    EXPECT_LT(v2_1, v2_3);
    EXPECT_LT(v2_2, v2_3);
    EXPECT_FALSE(v2_1 < v2_2);

    EXPECT_GT(v2_3, v2_0); // greater-than (all elements must be greater-than)
    EXPECT_GT(v2_3, v2_1);
    EXPECT_GT(v2_3, v2_2);
    EXPECT_FALSE(v2_2 > v2_1);

    // equality
    EXPECT_EQ(v2_3, v2_4);
    EXPECT_EQ(v3_3, v3_4);

    // Use the std::array for lexicographical order
    using Arr3 = typename Vec3::ArrayType;

    v3_1[1] = 4;
    const Arr3 &a3_0 = v3_0.getArray();
    const Arr3 &a3_1 = v3_1.getArray();
    const Arr3 &a3_2 = v3_2.getArray();
    const Arr3 &a3_3 = v3_3.getArray();
    const Arr3 &a3_4 = v3_4.getArray();

    // std::array
    EXPECT_NE(a3_0, a3_1);
    EXPECT_NE(a3_1, a3_2);
    EXPECT_NE(a3_2, a3_3);

    EXPECT_LT(a3_0, a3_1);
    EXPECT_LT(a3_1, a3_2);
    EXPECT_LT(a3_2, a3_3);

    EXPECT_LE(a3_0, a3_1);
    EXPECT_LE(a3_1, a3_2);
    EXPECT_LE(a3_2, a3_3);
    EXPECT_LE(a3_3, a3_4);

    EXPECT_GT(a3_1, a3_0);
    EXPECT_GT(a3_2, a3_1);
    EXPECT_GT(a3_3, a3_2);

    EXPECT_GE(a3_1, a3_0);
    EXPECT_GE(a3_2, a3_1);
    EXPECT_GE(a3_3, a3_2);
    EXPECT_GE(a3_3, a3_4);

    // Cubism::Core::Vector using lexicographic compare similar to
    // std::lexicographical_compare
    EXPECT_TRUE(v3_0.lexLT(v3_1));
    EXPECT_TRUE(v3_1.lexLT(v3_2));
    EXPECT_TRUE(v3_2.lexLT(v3_3));

    EXPECT_TRUE(v3_0.lexLE(v3_1));
    EXPECT_TRUE(v3_1.lexLE(v3_2));
    EXPECT_TRUE(v3_2.lexLE(v3_3));
    EXPECT_TRUE(v3_3.lexLE(v3_4));

    EXPECT_TRUE(v3_1.lexGT(v3_0));
    EXPECT_TRUE(v3_2.lexGT(v3_1));
    EXPECT_TRUE(v3_3.lexGT(v3_2));

    EXPECT_TRUE(v3_1.lexGE(v3_0));
    EXPECT_TRUE(v3_2.lexGE(v3_1));
    EXPECT_TRUE(v3_3.lexGE(v3_2));
    EXPECT_TRUE(v3_4.lexGE(v3_3));
}

// Cast
TEST(Vector, Cast)
{
    using Vec = V_f<3>;
    using Arr = typename Vec::ArrayType;
    using DataType = typename Vec::DataType;

    Vec v0;

    // cast to underlying array
    Arr a0 = static_cast<Arr>(v0);
    EXPECT_EQ(a0, v0.getArray());

    // cast to pointer to underlying data
    DataType *p0 = static_cast<DataType *>(v0);
    EXPECT_EQ(p0, v0.data());
}

// Arithmetic operators
TEST(Vector, Arithmetic)
{
    // it should be clear that some operators may yield undefined behavior for
    // unsigned types.
    using Vec = V_i<3>; // test with signed integers
    using DataType = typename Vec::DataType;

    const Vec v0;
    const Vec v1{1, 1, 1};
    const Vec v2{2, 2, 2};

    Vec aux0(v2);
    aux0 -= v1;
    ASSERT_EQ(aux0, v1);
    aux0 += v1;
    ASSERT_EQ(aux0, v2);
    aux0 /= v2;
    ASSERT_EQ(aux0, v1);
    aux0 *= v2;
    EXPECT_EQ(aux0, v2);
    aux0 = v1 + v1;
    EXPECT_EQ(aux0, v2);
    aux0 = v2 - v1;
    EXPECT_EQ(aux0, v1);
    aux0 = v2 / v2;
    EXPECT_EQ(aux0, v1);
    aux0 = v1 * v2;
    EXPECT_EQ(aux0, v2);

    const DataType fac1 = 1;
    const DataType fac2 = 2;
    aux0 = v1;
    aux0 += fac1;
    ASSERT_EQ(aux0, v2);
    aux0 -= fac1;
    ASSERT_EQ(aux0, v1);
    aux0 *= fac2;
    ASSERT_EQ(aux0, v2);
    aux0 /= fac2;
    EXPECT_EQ(aux0, v1);
    aux0 = v1 + fac1;
    EXPECT_EQ(aux0, v2);
    aux0 = v2 - fac1;
    EXPECT_EQ(aux0, v1);
    aux0 = v1 * fac2;
    EXPECT_EQ(aux0, v2);
    aux0 = v2 / fac2;
    EXPECT_EQ(aux0, v1);
    aux0 = fac1 + v1;
    EXPECT_EQ(aux0, v2);
    aux0 = fac2 - v1;
    EXPECT_EQ(aux0, v1);
    aux0 = fac2 * v1;
    EXPECT_EQ(aux0, v2);
    aux0 = fac2 / v2; // element wise scaled reciprocal: v2[i] = fac2 / v2[i]
    EXPECT_EQ(aux0, v1);

    aux0 = -1;
    EXPECT_EQ(aux0, -v1);
}

// Common vector operations (float and integral types)
TEST(Vector, CommonVecOp)
{
    using Vec = V_i<3>; // test with signed integers
    using DataType = typename Vec::DataType;

    const DataType nul = 0;
    const DataType one = 1;
    const DataType two = 2;
    const DataType tre = 3;
    Vec v0(two);
    Vec v1(tre);
    const DataType N = v0.size();

    EXPECT_EQ(v0.normL1(), two * N);
    EXPECT_EQ(v0.normLinf(), two);
    EXPECT_EQ(v0.normsq(), two * two * N);
    EXPECT_EQ(v0.dot(v1), two * tre * N);
    EXPECT_EQ(v0.distsq(v1), N); // works with unsigned integral for this test
    EXPECT_EQ(v0.sum(), two * N);
    EXPECT_EQ(v0.prod(), std::pow(two, N));
    EXPECT_EQ(v0.min(), two);
    EXPECT_EQ(v0.max(), two);
    EXPECT_EQ(v0.argmin(), static_cast<size_t>(N - one));
    EXPECT_EQ(v0.argmax(), static_cast<size_t>(N - one));

    v0 = Vec{1, 0};
    v1 = Vec{0, 1};
    EXPECT_EQ(v0.dot(v1), nul);

    for (size_t i = 0; i < v0.size(); ++i) {
        v0[i] = i + 1;
    }
    EXPECT_EQ(v0.min(), one);
    EXPECT_EQ(v0.max(), N);
    EXPECT_EQ(v0.argmin(), static_cast<size_t>(0));
    EXPECT_EQ(v0.argmax(), static_cast<size_t>(N - one));

    v0 = two;
    EXPECT_EQ(sumVector(v0.abs()), two * N);
    v0 = -two;
    EXPECT_EQ(sumVector(v0.abs()), two * N);
}

// Common vector operations with floating point return type
TEST(Vector, CommonVecOpReal)
{
    using Vec = V_d<3>;
    using DataType = typename Vec::DataType;

    const DataType one = 1;
    const DataType two = 2;
    const DataType tre = 3;
    Vec v0(two);
    Vec v1(tre);
    Vec vOne(one);

    EXPECT_EQ(v0.normL2(), std::sqrt(v0.dot(v0)));
    EXPECT_EQ(v0.norm(), std::sqrt(v0.dot(v0))); // alias for normL2
    EXPECT_EQ(v0.dist(v1), vOne.norm());
}

// Cross-product for 3D and 2D Cubism::Core::Vector
TEST(Vector, CrossProduct)
{
    using Vec2 = V_f<2>;
    using Vec3 = V_f<3>;
    using DataType = typename Vec3::DataType;

    // 3D Cubism::Core::Vector
    Vec3 v0{1, 0};
    Vec3 v1{0, 1};
    Vec3 v2 = v0.cross(v1);
    EXPECT_EQ(v2, Vec3({0, 0, 1}));
    v0[1] = 1;
    v1[0] = 1;
    v2 = v0.cross(v1);
    EXPECT_EQ(v2, Vec3({0}));

    // 2D Cubism::Core::Vector
    Vec2 v3{1, 0};
    Vec2 v4{0, 1};
    EXPECT_EQ(v3.getCrossThird(v4), static_cast<DataType>(1));
}

// Iterators
TEST(Vector, Iterator)
{
    using Vec = V_d<16>;
    using DataType = typename Vec::DataType;

    Vec v0;

    // forward
    for (DataType &v : v0) {
        v = 1;
    }
    EXPECT_EQ(sumVector(v0), v0.size());

    DataType sum = 0;
    for (const DataType &v : v0) {
        // v += 1; // compile-time error
        sum += v;
    }
    EXPECT_EQ(sum, v0.size());

    // reverse
    for (auto v = v0.rbegin(); v != v0.rend(); ++v) {
        *v += 1;
    }
    EXPECT_EQ(sumVector(v0), 2 * v0.size());

    sum = 0;
    for (auto v = v0.crbegin(); v != v0.crend(); ++v) {
        // *v += 1; // compile-time error
        sum += *v;
    }
    EXPECT_EQ(sum, 2 * v0.size());
}

} // namespace
