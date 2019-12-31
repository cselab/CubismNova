// File       : Field.h
// Created    : Mon Dec 30 2019 05:52:21 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data field
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef FIELD_H_7NZ0QFMC
#define FIELD_H_7NZ0QFMC

#include "Alloc/AlignedBlockAllocator.h"
#include "Block/Data.h"
#include "Core/Common.h"

#include <array>
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

template <typename TBlockData>
class Field : public TBlockData
{
protected:
    using TBlockData::block_;
    using TBlockData::range_;

    template <typename T>
    class IteratorBase
    {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T *;
        using reference = T &;

        IteratorBase(pointer ptr) : ptr_(ptr) {}
        IteratorBase(const IteratorBase &c) = default;
        ~IteratorBase() = default;

        IteratorBase &operator=(const IteratorBase &c) = default;
        IteratorBase &operator=(pointer ptr)
        {
            ptr_ = ptr;
            return *this;
        }

        bool operator==(const IteratorBase &rhs) const
        {
            return ptr_ == rhs.ptr_;
        }
        bool operator!=(const IteratorBase &rhs) const
        {
            return ptr_ != rhs.ptr_;
        }
        IteratorBase &operator+=(const difference_type rhs)
        {
            ptr_ += rhs;
            return *this;
        }
        IteratorBase &operator-=(const difference_type rhs)
        {
            ptr_ -= rhs;
            return *this;
        }
        IteratorBase &operator++()
        {
            ++ptr_;
            return *this;
        }
        IteratorBase &operator--()
        {
            --ptr_;
            return *this;
        }
        IteratorBase operator++(int)
        {
            IteratorBase tmp(*this);
            ++ptr_;
            return tmp;
        }
        IteratorBase operator--(int)
        {
            IteratorBase tmp(*this);
            --ptr_;
            return tmp;
        }
        IteratorBase operator+(const difference_type rhs)
        {
            IteratorBase tmp(*this);
            return (tmp += rhs);
        }
        IteratorBase operator-(const difference_type rhs)
        {
            IteratorBase tmp(*this);
            return (tmp -= rhs);
        }
        difference_type operator-(const IteratorBase &rhs) const
        {
            return std::distance(ptr_, rhs.ptr_);
        }
        reference operator*() { return *ptr_; }
        pointer operator->() { return ptr_; }

    private:
        pointer ptr_;
    };

public:
    using BaseType = TBlockData;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

    explicit Field(const IndexRangeType &r,
                   const size_t rank = 0,
                   const size_t comp = 0)
        : BaseType(r), nelements_(range_.size()), rank_(rank), comp_(comp)
    {
    }

    Field(const Field &f,
          const typename BaseType::MemoryOwner o,
          const size_t rank = 0,
          const size_t comp = 0)
        : BaseType(f, o), nelements_(range_.size()), rank_(rank), comp_(comp)
    {
    }

    Field(const IndexRangeType &r,
          DataType *ptr,
          const size_t bytes,
          const size_t rank = 0,
          const size_t comp = 0)
        : BaseType(r, ptr, bytes), nelements_(range_.size()), rank_(rank),
          comp_(comp)
    {
    }

    // The constructor
    //
    // Field(const IndexRangeType &r,
    //       const std::vector<DataType *> &ptr_list,
    //       const size_t bytes);
    //
    // is not implemented here on purpose.  It does not make sense for a scalar
    // field and must yield a compile time error when called.

    Field() = delete;
    Field(const Field &c) = default;
    Field(Field &&c) noexcept = default;
    Field &operator=(const Field &c) = default;
    Field &operator=(Field &&c) = default;
    ~Field() = default;

    using iterator = IteratorBase<DataType>;
    using const_iterator = IteratorBase<const DataType>;
    iterator begin() { return iterator(block_); }
    iterator end() { return iterator(block_ + range_.size()); }
    const_iterator cbegin() const { return const_iterator(block_); }
    const_iterator cend() const { return const_iterator(block_ + nelements_); }

    bool isScalar() const { return (0 == rank_); }
    size_t getRank() const { return rank_; }
    size_t getComp() const { return comp_; }

// TODO: [fabianw@mavt.ethz.ch; 2019-12-31]
#if 0
    Field operator-() const
    {
        Field f(*this);
        // fieldMul(f.block_, -1, f.block_, nelements_);
        for (size_t i = 0; i < nelements_; ++i) {
            f[i] = -f[i];
        }
        return f;
    }

    Field &operator+=(const Field &rhs)
    {
        fieldAdd(block_, rhs.block_, block_, nelements_);
        return *this;
    }
    Field &operator-=(const Field &rhs)
    {
        fieldSub(block_, rhs.block_, block_, nelements_);
        return *this;
    }
    Field &operator*=(const Field &rhs)
    {
        // element-wise multiplication
        fieldMul(block_, rhs.block_, block_, nelements_);
        return *this;
    }
    Field &operator/=(const Field &rhs)
    {
        // element-wise division
        fieldDiv(block_, rhs.block_, block_, nelements_);
        return *this;
    }

    // The following use RVO on non-const lhs
    friend Field operator+(Field lhs, const Field &rhs) { return (lhs += rhs); }

    friend Field operator-(Field lhs, const Field &rhs) { return (lhs -= rhs); }

    friend Field operator*(Field lhs, const Field &rhs) { return (lhs *= rhs); }

    friend Field operator/(Field lhs, const Field &rhs) { return (lhs /= rhs); }

    // Scalar rhs
    Field &operator+=(const DataType rhs)
    {
        fieldAdd(block_, rhs, block_, nelements_);
        return *this;
    }

    Field &operator-=(const DataType rhs)
    {
        fieldSub(block_, rhs, block_, nelements_);
        return *this;
    }

    Field &operator*=(const DataType rhs)
    {
        fieldMul(block_, rhs, block_, nelements_);
        return *this;
    }

    Field &operator/=(const DataType rhs)
    {
        assert(rhs > 0 || rhs < 0);
        fieldDiv(block_, rhs, block_, nelements_);
        return *this;
    }

    friend Field operator+(Field lhs, const DataType rhs)
    {
        return (lhs += rhs);
    }

    friend Field operator-(Field lhs, const DataType rhs)
    {
        return (lhs -= rhs);
    }

    friend Field operator*(Field lhs, const DataType rhs)
    {
        return (lhs *= rhs);
    }

    friend Field operator/(Field lhs, const DataType rhs)
    {
        return (lhs /= rhs);
    }

    friend Field operator+(const DataType lhs, Field rhs)
    {
        return (rhs += lhs);
    }

    friend Field operator-(const DataType lhs, Field rhs)
    {
        return -(rhs -= lhs);
    }

    friend Field operator*(const DataType lhs, Field rhs)
    {
        return (rhs *= lhs);
    }

    friend Field operator/(const DataType lhs, Field rhs)
    {
        // reciprocal multiplied by lhs
        fieldRcp(rhs.block_, lhs, rhs.block_, nelements_);
        return rhs;
    }
#endif

private:
    const size_t nelements_; // number of elements carried in block data
    const size_t rank_; // field rank -- rank_=0 (scalar), rank_=1 (vector), ...
    const size_t comp_; // field component in rank_ dimension
};

template <typename T, template <typename> class Alloc = AlignedBlockAllocator>
using CellField = Field<Data<T, DataMapping::Cell, CUBISM_DIMENSION, Alloc<T>>>;

template <typename T, template <typename> class Alloc = AlignedBlockAllocator>
using NodeField = Field<Data<T, DataMapping::Node, CUBISM_DIMENSION, Alloc<T>>>;

template <typename T, template <typename> class Alloc = AlignedBlockAllocator>
using FaceField = Field<Data<T, DataMapping::Face, CUBISM_DIMENSION, Alloc<T>>>;

#define FIELD_CONTAINER_OP_FIELD(OP)                                           \
    do {                                                                       \
        for (size_t i = 0; i < components_.size(); ++i) {                      \
            BaseType *LHS = components_[i];                                    \
            BaseType *RHS = rhs.components_[i];                                \
            if (LHS && RHS) {                                                  \
                *LHS OP *RHS;                                                  \
            }                                                                  \
        }                                                                      \
    } while (0)

#define FIELD_CONTAINER_OP_SCALAR(OP)                                          \
    do {                                                                       \
        for (size_t i = 0; i < components_.size(); ++i) {                      \
            BaseType *LHS = components_[i];                                    \
            if (LHS) {                                                         \
                *LHS OP rhs;                                                   \
            }                                                                  \
        }                                                                      \
    } while (0)

#define FIELD_CONTAINER_RCP()                                                  \
    do {                                                                       \
        for (size_t i = 0; i < rhs.components_.size(); ++i) {                  \
            BaseType *RHS = rhs.components_[i];                                \
            if (RHS) {                                                         \
                lhs / *RHS;                                                    \
            }                                                                  \
        }                                                                      \
    } while (0)

template <typename TField>
class FieldContainer
{
protected:
    using MemoryOwner = typename TField::BaseType::MemoryOwner;

public:
    using BaseType = TField;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

    explicit FieldContainer(const size_t n = 0) : components_(n, nullptr) {}

    FieldContainer(const size_t n,
                   const IndexRangeType &r,
                   const size_t rank = 0)
        : components_(n, nullptr)
    {
        for (size_t i = 0; i < n; ++i) {
            components_[i] = new BaseType(r, rank, i);
        }
    }

    FieldContainer(const std::vector<BaseType *> &vf)
        : components_(vf.size(), nullptr)
    {
        for (size_t i = 0; i < vf.size(); ++i) {
            const BaseType *comp = vf[i];
            if (comp) {
                components_[i] = new BaseType(*comp);
            }
        }
    }

    FieldContainer(const FieldContainer &fc,
                   const MemoryOwner owner,
                   const size_t rank = 0)
        : components_(fc.size(), nullptr)
    {
        for (size_t i = 0; i < fc.size(); ++i) {
            const BaseType *comp = fc.components_[i];
            if (comp) {
                BaseType *copy_comp;
                if (owner == MemoryOwner::Yes) {
                    copy_comp = new BaseType(*comp, owner, rank, i);
                } else {
                    const IndexRangeType r = comp->getIndexRange();
                    const size_t bytes = comp->getBlockBytes();
                    copy_comp = new BaseType(r, comp, bytes, rank, i);
                }
                components_[i] = copy_comp;
            }
        }
    }

    // The constructor
    //
    // FieldContainer(const IndexRangeType &r,
    //                const DataType *ptr,
    //                const size_t bytes);
    //
    // is not implemented here on purpose.  It does not make sense for a
    // container.  A compile time error will occur when this constructor is
    // called here.

    FieldContainer(const IndexRangeType &r,
                   const std::vector<DataType *> &ptr_list,
                   const size_t bytes,
                   const size_t rank = 0)
        : components_(ptr_list.size(), nullptr)
    {
        for (size_t i = 0; i < ptr_list.size(); ++i) {
            const BaseType *comp = ptr_list[i];
            if (comp) {
                components_[i] = new BaseType(r, comp, bytes, rank, i);
            }
        }
    }

    FieldContainer(const FieldContainer &c) : components_(c.size(), nullptr)
    {
        for (size_t i = 0; i < c.size(); ++i) {
            const BaseType *comp = c.components_[i];
            if (comp) {
                components_[i] = new BaseType(*comp);
            }
        }
    }

    FieldContainer(FieldContainer &&c) noexcept
        : components_(std::move(c.components_))
    {
        c.components_.clear();
    }

    virtual ~FieldContainer() { dispose_(); }

    FieldContainer &operator=(const FieldContainer &rhs)
    {
        if (this != &rhs) {
            dispose_();
            components_.resize(rhs.size());
            for (size_t i = 0; i < rhs.size(); ++i) {
                components_[i] = nullptr;
                BaseType *comp = rhs.components_[i];
                if (comp) {
                    components_[i] = new BaseType(*comp);
                }
            }
        }
        return *this;
    }

    FieldContainer &operator=(FieldContainer &&rhs)
    {
        if (this != &rhs) {
            dispose_();
            components_ = std::move(rhs.components_);
            rhs.components_.clear();
        }
        return *this;
    }

    size_t size() const { return components_.size(); }

    BaseType &operator[](const size_t i) { return getComp_(i); }
    const BaseType &operator[](const size_t i) const { return getComp_(i); }

    template <typename T>
    BaseType &operator[](const T &t)
    {
        return getComp_(static_cast<size_t>(t));
    }
    template <typename T>
    const BaseType &operator[](const T &t) const
    {
        return getComp_(static_cast<size_t>(t));
    }

#if 0
    // arithmetic
    FieldContainer operator-() const
    {
        FieldContainer fc(*this);
        for (size_t i = 0; i < fc.components_.size(); ++i) {
            BaseType *comp = fc.components_[i];
            if (comp) {
                // FIXME: [fabianw@mavt.ethz.ch; 2019-12-31] this is
                // inefficient
                *comp = -(*comp);
            }
        }
        return fc;
    }

    FieldContainer &operator+=(const FieldContainer &rhs)
    {
        FIELD_CONTAINER_OP_FIELD(+=);
        return *this;
    }
    FieldContainer &operator-=(const FieldContainer &rhs)
    {
        FIELD_CONTAINER_OP_FIELD(-=);
        return *this;
    }
    FieldContainer &operator*=(const FieldContainer &rhs)
    {
        FIELD_CONTAINER_OP_FIELD(*=);
        return *this;
    }
    FieldContainer &operator/=(const FieldContainer &rhs)
    {
        FIELD_CONTAINER_OP_FIELD(/=);
        return *this;
    }

    // The following use RVO on non-const lhs
    friend FieldContainer operator+(FieldContainer lhs,
                                    const FieldContainer &rhs)
    {
        return (lhs += rhs);
    }

    friend FieldContainer operator-(FieldContainer lhs,
                                    const FieldContainer &rhs)
    {
        return (lhs -= rhs);
    }

    friend FieldContainer operator*(FieldContainer lhs,
                                    const FieldContainer &rhs)
    {
        return (lhs *= rhs);
    }

    friend FieldContainer operator/(FieldContainer lhs,
                                    const FieldContainer &rhs)
    {
        return (lhs /= rhs);
    }

    // Scalar rhs
    FieldContainer &operator+=(const DataType rhs)
    {
        FIELD_CONTAINER_OP_SCALAR(+=);
        return *this;
    }

    FieldContainer &operator-=(const DataType rhs)
    {
        FIELD_CONTAINER_OP_SCALAR(-=);
        return *this;
    }

    FieldContainer &operator*=(const DataType rhs)
    {
        FIELD_CONTAINER_OP_SCALAR(*=);
        return *this;
    }

    FieldContainer &operator/=(const DataType rhs)
    {
        assert(rhs > 0 || rhs < 0);
        FIELD_CONTAINER_OP_SCALAR(/=);
        return *this;
    }

    friend FieldContainer operator+(FieldContainer lhs, const DataType rhs)
    {
        return (lhs += rhs);
    }

    friend FieldContainer operator-(FieldContainer lhs, const DataType rhs)
    {
        return (lhs -= rhs);
    }

    friend FieldContainer operator*(FieldContainer lhs, const DataType rhs)
    {
        return (lhs *= rhs);
    }

    friend FieldContainer operator/(FieldContainer lhs, const DataType rhs)
    {
        return (lhs /= rhs);
    }

    friend FieldContainer operator+(const DataType lhs, FieldContainer rhs)
    {
        return (rhs += lhs);
    }

    friend FieldContainer operator-(const DataType lhs, FieldContainer rhs)
    {
        return -(rhs -= lhs);
    }

    friend FieldContainer operator*(const DataType lhs, FieldContainer rhs)
    {
        return (rhs *= lhs);
    }

    friend FieldContainer operator/(const DataType lhs, FieldContainer rhs)
    {
        // reciprocal multiplied by lhs
        FIELD_CONTAINER_RCP();
        return rhs;
    }
#endif

protected:
    std::vector<BaseType *> components_;

private:
    BaseType &getComp_(const size_t i)
    {
        assert(i < components_.size());
        BaseType *comp = components_[i];
        if (!comp) {
            throw std::runtime_error("FieldContainer: Component " +
                                     std::to_string(i) +
                                     " was not assigned (nullptr)");
        }
        return *comp;
    }

    const BaseType &getComp_(const size_t i) const
    {
        assert(i < components_.size());
        const BaseType *comp = components_[i];
        if (!comp) {
            throw std::runtime_error("FieldContainer: Component " +
                                     std::to_string(i) +
                                     " was not assigned (nullptr)");
        }
        return *comp;
    }

    void dispose_()
    {
        for (size_t i = 0; i < components_.size(); ++i) {
            BaseType *comp = components_[i];
            if (comp) {
                delete comp;
                components_[i] = nullptr;
            }
        }
    }
};

#undef FIELD_CONTAINER_OP_FIELD
#undef FIELD_CONTAINER_OP_SCALAR
#undef FIELD_CONTAINER_RCP

template <typename T>
class FaceFieldContainer : public FieldContainer<FaceField<T>>
{
public:
    using FaceFieldType = FaceField<T>;
    using BaseType = FieldContainer<FaceFieldType>;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

    FaceFieldContainer() {}
    ~FaceFieldContainer() {}

private:
};

template <typename TField, size_t RANK>
class TensorField : public FieldContainer<TField>
{
    template <size_t B, size_t E>
    struct Power {
        static constexpr size_t value = B * Power<B, E - 1>::value;
    };
    template <size_t B>
    struct Power<B, 0> {
        static constexpr size_t value = 1;
    };

    using MemoryOwner = typename TField::BaseType::MemoryOwner;

public:
    using BaseType = TField;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

    static constexpr size_t Rank = RANK;
    static constexpr size_t NComponents =
        Power<IndexRangeType::Dim, RANK>::value;

    explicit TensorField(const IndexRangeType &r)
        : FieldContainer<TField>(NComponents, r, Rank)
    {
    }

    TensorField(const TensorField &fc, const MemoryOwner owner)
        : FieldContainer<TField>(fc, owner, Rank)
    {
    }

    TensorField(const IndexRangeType &r,
                const std::vector<DataType *> &ptr_list,
                const size_t bytes)
        : FieldContainer<TField>(r, ptr_list, bytes, Rank)
    {
    }

    TensorField(const TensorField &c) = default;
    TensorField(TensorField &&c) noexcept = default;
    TensorField &operator=(const TensorField &c) = default;
    TensorField &operator=(TensorField &&c) = default;
    ~TensorField() = default;
};

template <typename TField, size_t RANK>
constexpr size_t TensorField<TField, RANK>::Rank;

template <typename TField, size_t RANK>
constexpr size_t TensorField<TField, RANK>::NComponents;

template <typename TField>
using VectorField = TensorField<TField, 1>;

/// @brief Field view type (never owns block memory)
///
/// @tparam TField (underlying field type)
template <typename TField>
class FieldView : public TField
{
    using MemoryOwner = typename TField::BaseType::MemoryOwner;

public:
    using BaseType = TField;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

    FieldView(const BaseType &f) : BaseType(f, MemoryOwner::No) {}

    FieldView(const IndexRangeType &r, DataType *ptr, const size_t bytes)
        : BaseType(r, ptr, bytes)
    {
    }

    FieldView(const IndexRangeType &r,
              const std::vector<DataType *> &ptr_list,
              const size_t bytes)
        : BaseType(r, ptr_list, bytes)
    {
    }

    FieldView() = delete;
    FieldView(const FieldView &c) = default;
    FieldView(FieldView &&c) noexcept = default;
    FieldView &operator=(const FieldView &c) = default;
    FieldView &operator=(FieldView &&c) = default;

    ~FieldView() = default;
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* FIELD_H_7NZ0QFMC */
