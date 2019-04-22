// File       : Vector.h
// Created    : Mon Apr 22 2019 11:36:26 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Genereic vector class with support for common operations
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef VECTOR_H_8YBMEXHP
#define VECTOR_H_8YBMEXHP

#include "Core/Common.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <initializer_list>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

NAMESPACE_BEGIN(Cubism)

/// @brief Vector class with support for common operations (wraps around
///        std::array)
///
/// @tparam T   (underlying data type, must follow the rules of aggreagte
///             initialization)
/// @tparam DIM (vector dimension, should be low-dimensional when used for
///             automatic variables on the stack)
template <typename T, size_t DIM>
class Vector
{
public:
    using DataType = T;
    using HashType = std::array<unsigned char, DIM * sizeof(T)>;

    static constexpr size_t Dim = DIM;

    /// @brief Compute i-th vector component from byte hash h
    ///
    /// @param i
    /// @param h
    ///
    /// @return
    static DataType getComponent(const size_t i, const HashType &h)
    {
        assert(i < DIM);
        const DataType *c = reinterpret_cast<const DataType *>(h.data());
        return c[i];
    }

    /// @brief Compute hex-string representation of byte hash
    ///
    /// @param h
    ///
    /// @return
    static std::string getHashString(const HashType &h)
    {
        std::ostringstream os;
        os << "0x";
        for (auto rit = h.rbegin(); rit < h.rend(); ++rit) {
            os << std::hex << std::setfill('0') << std::setw(2)
               << static_cast<int>(*rit);
        }
        return os.str();
    }

    /// @brief Default constructor
    Vector() : hash_({0}), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
    }

    /// @brief Default copy constructor
    Vector(const Vector &c)
        : hash_(c.hash_), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
    }

    /// @brief Default move constructor
    Vector(Vector &&c) noexcept
        : hash_(std::move(c.hash_)),
          data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
        c.data_ = nullptr;
    }

    /// @brief Construct vector from byte hash
    Vector(const HashType &h)
        : hash_(h), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
    }

    /// @brief Constructor for arbitrary Vector types
    ///
    /// @tparam U
    /// @tparam DIMU
    /// @param c
    template <typename U, size_t DIMU>
    Vector(const Vector<U, DIMU> &c)
        : hash_({0}), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
        copy_from_address_(c.data(), c.size());
    }

    /// @brief Constructor for any scalar type U.  The type U must be castable
    ///        to DataType.
    ///
    /// @tparam U
    /// @param scalar
    template <typename U>
    explicit Vector(const U scalar)
        : hash_({0}), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
        const DataType v = static_cast<DataType>(scalar);
        std::fill(data_, data_ + DIM, v);
    }

    /// @brief Constructor for list initialization
    ///
    /// @tparam U
    /// @param ilist
    template <typename U>
    Vector(std::initializer_list<U> ilist)
        : hash_({0}), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
        copy_from_address_(ilist.begin(), ilist.size());
    }

    /// @brief Constructor to initialize data from std::vector
    ///
    /// @tparam U
    /// @param vec
    template <typename U>
    Vector(const std::vector<U> &vec)
        : hash_({0}), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
        copy_from_address_(vec.data(), vec.size());
    }

    /// @brief Constructor to initialize data from std::array
    ///
    /// @tparam U
    /// @tparam DIMU
    /// @param ary
    template <typename U, size_t DIMU>
    Vector(const std::array<U, DIMU> &ary)
        : hash_({0}), data_(reinterpret_cast<DataType *>(hash_.data()))
    {
        assert(static_cast<void *>(&hash_[0]) == static_cast<void *>(data_));
        copy_from_address_(ary.data(), ary.size());
    }

    /// @brief Default destructor
    virtual ~Vector() {}

    /// @brief Default assignment operator
    Vector &operator=(const Vector &c)
    {
        if (this != &c) {
            hash_ = c.hash_;
        }
        return *this;
    }

    /// @brief Default move assignment operator
    Vector &operator=(Vector &&c)
    {
        if (this != &c) {
            hash_ = std::move(c.hash_);
            data_ = reinterpret_cast<DataType *>(hash_.data());
            c.data_ = nullptr;
            assert(static_cast<void *>(&hash_[0]) ==
                   static_cast<void *>(data_));
        }
        return *this;
    }

    /// @brief Assigment operator for any other Vector type.  The data type U
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
            copy_from_address_(c.data(), c.size());
        }
        return *this;
    }

    /// @brief Swap this vector with other vector
    ///
    /// @param other
    void swap(Vector &other)
    {
        using std::swap;
        swap(hash_, other.hash_);
        // std::array does not swap the memory addresses.  the following two
        // lines are a safety measure nevertheless.
        data_ = reinterpret_cast<DataType *>(hash_.data());
        other.data_ = reinterpret_cast<DataType *>(other.hash_.data());
    }

    /// @brief Return byte hash
    const HashType &getHash() const { return hash_; }

    /// @brief Return size of vector
    size_t size() const { return DIM; }

    /// @brief Return raw data
    DataType *data() { return data_; }
    const DataType *data() const { return data_; }

    /// @brief Data access interface
    DataType &operator[](const size_t i)
    {
        assert(i < DIM);
        return data_[i];
    }

    /// @brief Data access interface
    const DataType &operator[](const size_t i) const
    {
        assert(i < DIM);
        return data_[i];
    }

    // Arithmetic operators

    // Common vector operations interface

private:
    HashType hash_;  // underlying data represented in bytes
    DataType *data_; // data represented as DataType

    template <typename U>
    void copy_from_address_(const U *src, size_t size_src)
    {
        // Scalar type U must be castable to DataType.  If size_src < DIM then
        // the remaining elements will be left uninitialized.
        if (size_src > DIM) {
            size_src = DIM;
        }
        DataType *dst = data_;
        for (size_t i = 0; i < size_src; ++i) {
            *dst++ = static_cast<DataType>(*src++);
        }
    }
};

template <typename T, size_t DIM>
void swap(Vector<T, DIM> &va, Vector<T, DIM> &vb)
{
    va.swap(vb);
}

NAMESPACE_END(Cubism)

#endif /* VECTOR_H_8YBMEXHP */
