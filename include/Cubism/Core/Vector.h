// File       : Vector.h
// Created    : Mon Apr 22 2019 11:36:26 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Generic vector class with support for common operations
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef VECTOR_H_8YBMEXHP
#define VECTOR_H_8YBMEXHP

#include "Core/Common.h"
#include "Core/Math.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

NAMESPACE_BEGIN(Cubism)

/// @brief Vector class with support for common operations (wraps around
///        std::array)
///
/// @tparam T   (underlying data type, must follow the rules of aggregate
///             initialization)
/// @tparam DIM (vector dimension, should be low-dimensional when used for
///             automatic variables on the stack)
template <typename T, size_t DIM>
class Vector
{
public:
    using ArrayType = std::array<T, DIM>;
    using DataType = T;

    static constexpr size_t Bytes = DIM * sizeof(T);
    static constexpr size_t Dim = DIM;

    using iterator = typename ArrayType::iterator;
    using const_iterator = typename ArrayType::const_iterator;
    using reverse_iterator = typename ArrayType::reverse_iterator;
    using const_reverse_iterator = typename ArrayType::const_reverse_iterator;

    template <size_t DIR>
    static Vector getUnitVector()
    {
        assert(DIR < Dim);
        Vector u;
        u[DIR] = 1;
        return u;
    }

    /// @brief Default constructor
    Vector() : array_({0}) {}

    /// @brief Default copy constructor
    Vector(const Vector &c) : array_(c.array_) {}

    /// @brief Default move constructor (linear complexity)
    Vector(Vector &&c) noexcept : array_(std::move(c.array_)) {}

    /// @brief Copy constructor to initialize data from std::array of same type
    Vector(const ArrayType &ary) : array_(ary) {}

    /// @brief Move constructor to initialize data from std::array of same type
    ///        (linear complexity)
    Vector(ArrayType &&ary) noexcept : array_(std::move(ary)) {}

    // Arbitrary type constructors:
    //
    /// @brief Constructor for arbitrary Vector types
    template <typename U, size_t DIMU>
    Vector(const Vector<U, DIMU> &c) : array_({0})
    {
        copyFromAddress_(c.data(), c.size());
    }

    /// @brief Constructor for any scalar type U.  The type U must be castable
    ///        to DataType.
    template <typename U>
    explicit Vector(const U scalar) : array_({0})
    {
        const DataType v = static_cast<DataType>(scalar);
        DataType *dst = this->data();
        std::fill(dst, dst + Dim, v);
    }

    /// @brief Constructor for arbitrary list initialization
    template <typename U>
    Vector(std::initializer_list<U> ilist) : array_({0})
    {
        copyFromAddress_(ilist.begin(), ilist.size());
    }

    /// @brief Constructor to initialize data from std::vector
    template <typename U>
    Vector(const std::vector<U> &vec) : array_({0})
    {
        copyFromAddress_(vec.data(), vec.size());
    }

    /// @brief Constructor to initialize data from arbitrary std::array
    template <typename U, size_t DIMU>
    Vector(const std::array<U, DIMU> &ary) : array_({0})
    {
        copyFromAddress_(ary.data(), ary.size());
    }

    /// @brief Constructor to initialize data from arbitrary pointer ptr.  The
    ///        number of elements n must be provided.
    template <typename U>
    explicit Vector(const U *ptr, size_t n) : array_({0})
    {
        copyFromAddress_(ptr, n);
    }

    /// @brief Default destructor
    virtual ~Vector() = default;

    /// @brief Default assignment operator
    Vector &operator=(const Vector &c)
    {
        if (this != &c) {
            array_ = c.array_;
        }
        return *this;
    }

    /// @brief Default move assignment operator (linear complexity)
    Vector &operator=(Vector &&c)
    {
        if (this != &c) {
            array_ = std::move(c.array_);
        }
        return *this;
    }

    /// @brief ArrayType assignment operator
    Vector &operator=(const ArrayType &ary)
    {
        if (&array_ != &ary) {
            array_ = ary;
        }
        return *this;
    }

    /// @brief ArrayType move assignment operator (linear complexity)
    Vector &operator=(ArrayType &&ary)
    {
        if (&array_ != &ary) {
            array_ = std::move(ary);
        }
        return *this;
    }

    /// @brief Assignment operator for scalar c
    Vector &operator=(const DataType c)
    {
        DataType *dst = this->data();
        std::fill(dst, dst + Dim, c);
        return *this;
    }

    /// @brief Assignment operator for any other Vector type.  The data type U
    ///        must be castable to DataType.  This operation is less efficient
    ///        than assigning vectors of the same type.
    ///
    /// @tparam U (data type must be castable to DataType)
    /// @tparam DIMU (DIMU may not necessarily be equal to DIM.  If DIMU < DIM,
    ///              the remaining data elements are left unchanged.)
    /// @param c
    ///
    /// @return (Vector<T, DIM>)
    template <typename U, size_t DIMU>
    Vector &operator=(const Vector<U, DIMU> &c)
    {
        if (static_cast<const void *>(this) != static_cast<const void *>(&c)) {
            copyFromAddress_(c.data(), c.size());
        }
        return *this;
    }

    /// @brief Swap this vector with other vector
    void swap(Vector &other) noexcept { array_.swap(other.array_); }

    /// @brief Return size of vector
    size_t size() const { return Dim; }

    /// @brief Return raw data
    DataType *data() { return array_.data(); }
    const DataType *data() const { return array_.data(); }

    /// @brief Return underlying std::array
    ArrayType &getArray() { return array_; }
    const ArrayType &getArray() const { return array_; }

    /// @brief Data access interface
    DataType &operator[](const size_t i) { return array_[i]; }

    /// @brief Data access interface
    const DataType &operator[](const size_t i) const { return array_[i]; }

    // Allowed casts:
    //
    /// @brief Cast vector to its underlying std::array type
    explicit operator ArrayType() const { return array_; }

    /// @brief Cast vector to pointer, pointing to the first element of its data
    explicit operator DataType *() { return data(); }

    // iterators:
    //
    iterator begin() noexcept { return array_.begin(); }
    iterator end() noexcept { return array_.end(); }
    reverse_iterator rbegin() noexcept { return array_.rbegin(); }
    reverse_iterator rend() noexcept { return array_.rend(); }

    const_iterator cbegin() const noexcept { return array_.cbegin(); }
    const_iterator cend() const noexcept { return array_.cend(); }
    const_reverse_iterator crbegin() const noexcept { return array_.crbegin(); }
    const_reverse_iterator crend() const noexcept { return array_.crend(); }

    // Comparison operators:
    //
    /// @brief Equality operator.  A Vector is equal to another iff they
    ///        are equal component-wise.
    bool operator==(const Vector &other) const
    {
        bool is_equal = true;
        for (size_t i = 0; i < Dim; ++i) {
            is_equal = is_equal && (!(array_[i] > other.array_[i]) &&
                                    !(array_[i] < other.array_[i]));
        }
        return is_equal;
    }
    bool operator!=(const Vector &other) const { return !(*this == other); }

    /// @brief Less than operator.  A Vector is smaller than another iff
    ///        all components are less than the corresponding components of
    ///        the other Vector.
    bool operator<(const Vector &other) const
    {
        bool is_less = true;
        for (size_t i = 0; i < Dim; ++i) {
            is_less = is_less && (array_[i] < other.array_[i]);
        }
        return is_less;
    }
    bool operator>(const Vector &other) const { return (other < *this); }

    /// @brief Less-equal than operator.  A Vector is smaller or equal
    ///        than another iff all components are less or equal than the
    ///        corresponding components of the other Vector.
    bool operator<=(const Vector &other) const
    {
        bool is_lesseq = true;
        for (size_t i = 0; i < Dim; ++i) {
            is_lesseq = is_lesseq && (array_[i] <= other.array_[i]);
        }
        return is_lesseq;
    }
    bool operator>=(const Vector &other) const { return (other <= *this); }

    // lexicographic compare
    bool lexLT(const Vector &other) const { return array_ < other.array_; }
    bool lexLE(const Vector &other) const { return array_ <= other.array_; }
    bool lexGT(const Vector &other) const { return array_ > other.array_; }
    bool lexGE(const Vector &other) const { return array_ >= other.array_; }

    // Unitary operators:
    //
    Vector operator-() const
    {
        Vector v(*this);
        for (size_t i = 0; i < Dim; ++i) {
            v[i] = -v[i];
        }
        return v;
    }

    // Arithmetic operators
    //
    // Vector rhs:
    Vector &operator+=(const Vector &rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] += rhs[i];
        }
        return *this;
    }

    Vector &operator-=(const Vector &rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] -= rhs[i];
        }
        return *this;
    }

    Vector &operator*=(const Vector &rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] *= rhs[i];
        }
        return *this;
    }

    Vector &operator/=(const Vector &rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            assert(rhs[i] > 0 || rhs[i] < 0);
            array_[i] /= rhs[i];
        }
        return *this;
    }

    // The following use RVO on non-const lhs
    friend Vector operator+(Vector lhs, const Vector &rhs)
    {
        return (lhs += rhs);
    }

    friend Vector operator-(Vector lhs, const Vector &rhs)
    {
        return (lhs -= rhs);
    }

    friend Vector operator*(Vector lhs, const Vector &rhs)
    {
        return (lhs *= rhs);
    }

    friend Vector operator/(Vector lhs, const Vector &rhs)
    {
        return (lhs /= rhs);
    }

    // Scalar rhs
    Vector &operator+=(const DataType rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] += rhs;
        }
        return *this;
    }

    Vector &operator-=(const DataType rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] -= rhs;
        }
        return *this;
    }

    Vector &operator*=(const DataType rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] *= rhs;
        }
        return *this;
    }

    Vector &operator/=(const DataType rhs)
    {
        assert(rhs > 0 || rhs < 0);
        for (size_t i = 0; i < Dim; ++i) {
            array_[i] /= rhs;
        }
        return *this;
    }

    friend Vector operator+(Vector lhs, const DataType rhs)
    {
        return (lhs += rhs);
    }

    friend Vector operator-(Vector lhs, const DataType rhs)
    {
        return (lhs -= rhs);
    }

    friend Vector operator*(Vector lhs, const DataType rhs)
    {
        return (lhs *= rhs);
    }

    friend Vector operator/(Vector lhs, const DataType rhs)
    {
        return (lhs /= rhs);
    }

    friend Vector operator+(const DataType lhs, Vector rhs)
    {
        return (rhs += lhs);
    }

    friend Vector operator-(const DataType lhs, Vector rhs)
    {
        // use single loop instead of two loops as in -(rhs -= lhs) [+ copy]
        for (size_t i = 0; i < Dim; ++i) {
            rhs[i] = lhs - rhs[i];
        }
        return rhs;
    }

    friend Vector operator*(const DataType lhs, Vector rhs)
    {
        return (rhs *= lhs);
    }

    friend Vector operator/(const DataType lhs, Vector rhs)
    {
        for (size_t i = 0; i < Dim; ++i) {
            assert(rhs[i] > 0 || rhs[i] < 0);
            rhs[i] = lhs / rhs[i];
        }
        return rhs;
    }

    // Common vector operations interface:
    //
    /// @brief Squared Euclidean vector norm
    DataType normsq() const { return sumProd_(*this); }

    /// @brief Euclidean vector norm (L2)
    DataType norm() const { return mySqrt(normsq()); }

    /// @brief Vector L2 norm (alias for norm())
    DataType normL2() const { return norm(); }

    /// @brief Vector L1 norm
    DataType normL1() const
    {
        DataType res = 0;
        for (size_t i = 0; i < Dim; ++i) {
            res += myAbs(array_[i]);
        }
        return res;
    }

    /// @brief Vector maximum norm
    DataType normLinf() const
    {
        DataType res = 0;
        for (size_t i = 0; i < Dim; ++i) {
            res = std::max(res, myAbs(array_[i]));
        }
        return res;
    }

    /// @brief Vector dot product
    DataType dot(const Vector &other) const { return sumProd_(other); }

    /// @brief Vector cross product
    Vector cross(const Vector &other) const
    {
        static_assert(3 == Dim, "Cross-product is not defined for Dim != 3");
        DataType v0 = array_[1] * other[2] - array_[2] * other[1];
        DataType v1 = array_[2] * other[0] - array_[0] * other[2];
        DataType v2 = array_[0] * other[1] - array_[1] * other[0];
        return Vector({v0, v1, v2});
    }

    /// @brief Third component of vector cross-product
    DataType getCrossThird(const Vector &other) const
    {
        static_assert(
            3 == Dim || 2 == Dim,
            "Third component of Cross-product requires Dim == 2 or Dim == 3");
        return array_[0] * other[1] - array_[1] * other[0];
    }

    /// @brief Squared distance between this and other vector
    DataType distsq(Vector other) const
    {
        other -= *this;
        return other.normsq();
    }

    /// @brief Distance between this and other vector
    DataType dist(Vector other) const
    {
        other -= *this;
        return other.norm();
    }

    /// @brief Sum of vector components
    DataType sum() const
    {
        DataType res = array_[0];
        for (size_t i = 1; i < Dim; ++i) {
            res += array_[i];
        }
        return res;
    }

    /// @brief Product of vector components
    DataType prod() const
    {
        DataType res = array_[0];
        for (size_t i = 1; i < Dim; ++i) {
            res *= array_[i];
        }
        return res;
    }

    /// @brief Minimum vector component
    DataType min() const
    {
        DataType res = array_[0];
        for (size_t i = 1; i < Dim; ++i) {
            res = std::min(res, array_[i]);
        }
        return res;
    }

    /// @brief Maximum vector component
    DataType max() const
    {
        DataType res = array_[0];
        for (size_t i = 1; i < Dim; ++i) {
            res = std::max(res, array_[i]);
        }
        return res;
    }

    /// @brief Index of minimum component
    size_t argmin() const
    {
        size_t res = 0;
        for (size_t i = 1; i < Dim; ++i) {
            res = (array_[res] < array_[i]) ? res : i;
        }
        return res;
    }

    /// @brief Index of maximum component
    size_t argmax() const
    {
        size_t res = 0;
        for (size_t i = 1; i < Dim; ++i) {
            res = (array_[res] > array_[i]) ? res : i;
        }
        return res;
    }

    /// @brief Copy of this vector with absolute value for all components
    Vector abs() const
    {
        Vector res(*this);
        for (size_t i = 0; i < Dim; ++i) {
            res[i] = myAbs(res[i]);
        }
        return res;
    }

private:
    ArrayType array_;

    /// @brief Copy from arbitrary source type.  The data type U must be
    ///        castable to DataType.  If size_src < Dim, then the remaining
    ///        elements will be left untouched.
    template <typename U>
    void copyFromAddress_(const U *src, size_t size_src)
    {
        if (size_src > Dim) {
            size_src = Dim;
        }
        for (size_t i = 0; i < size_src; ++i) {
            array_[i] = static_cast<DataType>(*src++);
        }
    }

    /// @brief Compute the sum of the component-wise products between this
    ///        vector and other vector.
    DataType sumProd_(const Vector &other) const
    {
        DataType res = 0;
        for (size_t i = 0; i < Dim; ++i) {
            res += array_[i] * other[i];
        }
        return res;
    }
};

/// @brief Non-STL swap function for vectors
template <typename T, size_t DIM>
void swap(Vector<T, DIM> &va, Vector<T, DIM> &vb) noexcept
{
    va.swap(vb);
}

/// @brief Input stream for Vector<T, DIM> types
template <typename T, size_t DIM>
std::istream &operator>>(std::istream &stream, Vector<T, DIM> &vec)
{
    assert(DIM > 0);
    for (size_t i = 0; i < DIM; ++i) {
        stream >> vec[i];
    }
    return stream;
}

/// @brief Output stream for Vector<T, DIM> types
template <typename T, size_t DIM>
std::ostream &operator<<(std::ostream &stream, const Vector<T, DIM> &vec)
{
    assert(DIM > 0);
    stream << "[" << vec[0];
    for (size_t i = 1; i < DIM; ++i) {
        stream << ", " << vec[i];
    }
    stream << "]";
    return stream;
}

NAMESPACE_END(Cubism)

#endif /* VECTOR_H_8YBMEXHP */
