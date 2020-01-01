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

struct FieldState {
    size_t rank;      // field rank -- rank=0 (scalar), rank=1 (vector), ...
    size_t comp;      // field component in rank dimension
};

// TODO: [fabianw@mavt.ethz.ch; 2020-01-01]
// class FieldUnit // for physical units of carried data

template <typename TBlockData, typename TState = FieldState>
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
    using FieldType = Field;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;
    using FieldStateType = TState;

    Field() = delete;

    explicit Field(const IndexRangeType &r,
                   const FieldStateType &fs = FieldStateType())
        : BaseType(r), state_(new FieldStateType())
    {
        *state_ = fs;
    }

    Field(const Field &f, const typename BaseType::MemoryOwner o)
        : BaseType(f, o)
    {
        if (this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(f);
    }

    Field(const IndexRangeType &r,
          DataType *ptr,
          const size_t bytes,
          FieldStateType *sptr)
        : BaseType(r, ptr, bytes), state_(sptr)
    {
    }

    Field(const Field &c) : BaseType(c)
    {
        if (this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(c);
    }

    Field(Field &&c) noexcept
        : BaseType(std::move(c)), state_(std::move(c.state_))
    {
        c.state_ = nullptr;
    }

    ~Field() { disposeState_(); }

    Field &operator=(const Field &c)
    {
        assert(range_.size() == c.range_.size());
        if (this != &c) {
            copyState_(c);
            BaseType::operator=(c);
        }
        return *this;
    };

    Field &operator=(Field &&c)
    {
        assert(range_.size() == c.range_.size());
        if (this != &c) {
            BaseType::operator=(std::move(c));
            disposeState_();
            state_ = std::move(c.state_);
            c.state_ = nullptr;
        }
        return *this;
    }

    using iterator = IteratorBase<DataType>;
    using const_iterator = IteratorBase<const DataType>;
    iterator begin() noexcept { return iterator(block_); }
    iterator end() noexcept { return iterator(block_ + range_.size()); }
    const_iterator cbegin() const noexcept { return const_iterator(block_); }
    const_iterator cend() const noexcept
    {
        return const_iterator(block_ + range_.size());
    }

    bool isScalar() const { return (0 == state_->rank); }
    size_t getRank() const { return state_->rank; }
    size_t getComp() const { return state_->comp; }
    FieldStateType &getState() { return *state_; }
    const FieldStateType &getState() const { return *state_; }

// TODO: [fabianw@mavt.ethz.ch; 2019-12-31]
#if 0
    Field operator-() const
    {
        Field f(*this);
        // fieldMul(f.block_, -1, f.block_, range_.size());
        for (size_t i = 0; i < range_.size(); ++i) {
            f[i] = -f[i];
        }
        return f;
    }

    Field &operator+=(const Field &rhs)
    {
        fieldAdd(block_, rhs.block_, block_, range_.size());
        return *this;
    }
    Field &operator-=(const Field &rhs)
    {
        fieldSub(block_, rhs.block_, block_, range_.size());
        return *this;
    }
    Field &operator*=(const Field &rhs)
    {
        // element-wise multiplication
        fieldMul(block_, rhs.block_, block_, range_.size());
        return *this;
    }
    Field &operator/=(const Field &rhs)
    {
        // element-wise division
        fieldDiv(block_, rhs.block_, block_, range_.size());
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
        fieldAdd(block_, rhs, block_, range_.size());
        return *this;
    }

    Field &operator-=(const DataType rhs)
    {
        fieldSub(block_, rhs, block_, range_.size());
        return *this;
    }

    Field &operator*=(const DataType rhs)
    {
        fieldMul(block_, rhs, block_, range_.size());
        return *this;
    }

    Field &operator/=(const DataType rhs)
    {
        assert(rhs > 0 || rhs < 0);
        fieldDiv(block_, rhs, block_, range_.size());
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
        fieldRcp(rhs.block_, lhs, rhs.block_, range_.size());
        return rhs;
    }
#endif

private:
    FieldStateType *state_; // Field state values/meta data

    void disposeState_()
    {
        if (this->isMemoryOwner()) {
            if (state_) {
                delete state_;
            }
        }
        state_ = nullptr;
    }

    void copyState_(const Field &c)
    {
        if (this->isMemoryOwner()) {
            *state_ = *c.state_;
        } else {
            state_ = c.state_;
        }
    }
};

template <typename T,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using CellField =
    Field<Data<T, DataMapping::Cell, CUBISM_DIMENSION, Alloc<T>>, State>;

template <typename T,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using NodeField =
    Field<Data<T, DataMapping::Node, CUBISM_DIMENSION, Alloc<T>>, State>;

template <typename T,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using FaceField =
    Field<Data<T, DataMapping::Face, CUBISM_DIMENSION, Alloc<T>>, State>;

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
                *RHS = lhs / *RHS;                                             \
            }                                                                  \
        }                                                                      \
    } while (0)

template <typename TField>
class FieldContainer
{
public:
    using BaseType = TField;
    using FieldType = typename BaseType::FieldType;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

private:
    using ContainerType = std::vector<BaseType *>;
    using FieldState = typename FieldType::FieldStateType;

protected:
    using MemoryOwner = typename TField::BaseType::MemoryOwner;

public:
    explicit FieldContainer(const size_t n = 0) : components_(n, nullptr) {}

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

    FieldContainer(const FieldContainer &fc, const MemoryOwner owner)
        : components_(fc.size(), nullptr)
    {
        for (size_t i = 0; i < fc.size(); ++i) {
            const BaseType *comp = fc.components_[i];
            if (comp) {
                components_[i] = new BaseType(*comp, owner);
            }
        }
    }

    FieldContainer(const size_t n,
                   const IndexRangeType &r,
                   const size_t rank = 0)
        : components_(n, nullptr)
    {
        FieldState fs;
        fs.rank = rank;
        for (size_t i = 0; i < n; ++i) {
            fs.comp = i;
            components_[i] = new FieldType(r, fs);
        }
    }

    FieldContainer(const std::vector<IndexRangeType> &range_list,
                   const std::vector<DataType *> &ptr_list,
                   const std::vector<size_t> &bytes_list,
                   const std::vector<FieldState *> &state_list)
        : components_(ptr_list.size(), nullptr)
    {
        assert(range_list.size() == ptr_list.size());
        assert(ptr_list.size() == bytes_list.size());
        for (size_t i = 0; i < ptr_list.size(); ++i) {
            DataType *data = ptr_list[i];
            assert(data != nullptr);
            components_[i] = new FieldType(
                range_list[i], data, bytes_list[i], state_list[i]);
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

    // WARNING: if you use an incomplete container (i.e. some components
    // are nullptr), the iterators will just return nullptr -- no further checks
    // are performed.
    using iterator = typename ContainerType::iterator;
    using const_iterator = typename ContainerType::const_iterator;
    using reverse_iterator = typename ContainerType::reverse_iterator;
    using const_reverse_iterator =
        typename ContainerType::const_reverse_iterator;
    iterator begin() noexcept { return components_.begin(); }
    iterator end() noexcept { return components_.end(); }
    reverse_iterator rbegin() noexcept { return components_.rbegin(); }
    reverse_iterator rend() noexcept { return components_.rend(); }
    const_iterator cbegin() const noexcept { return components_.cbegin(); }
    const_iterator cend() const noexcept { return components_.cend(); }
    const_reverse_iterator crbegin() const noexcept
    {
        return components_.crbegin();
    }
    const_reverse_iterator crend() const noexcept
    {
        return components_.crend();
    }

    size_t size() const { return components_.size(); }

    void copyData(const FieldContainer &rhs)
    {
        assert(components_.size() == rhs.components_.size());
        for (size_t i = 0; i < components_.size(); ++i) {
            BaseType *dst = components_[i];
            BaseType *src = rhs.components_[i];
            assert(dst && src);
            dst->copyData(*src);
        }
    }

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-01] push_back, remove (?)

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
    ContainerType components_;

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
class FaceFieldAll : public FieldContainer<FaceField<T>>
{
public:
    using BaseType = FieldContainer<FaceField<T>>;
    using FieldType = typename BaseType::FieldType;
    using FieldState = typename FieldType::FieldStateType;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

private:
    using MemoryOwner = typename FieldType::BaseType::MemoryOwner;

public:
    explicit FaceFieldAll(const IndexRangeType &cell_domain)
        : BaseType(IndexRangeType::Dim)
    {
        const MultiIndex cells = cell_domain.getExtent(); // number of cells
        FieldState fs;
        fs.rank = 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            // XXX: [fabianw@mavt.ethz.ch; 2020-01-01] Not most favorable for
            // vectorization but more intuitive.  Will require some
            // transposition in vectorized code.
            const IndexRangeType ri(cells + MultiIndex::getUnitVector(i));
            fs.comp = i;
            components_[i] = new FieldType(ri, fs);
            assert(components_[i] != nullptr);
        }
    }

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-01] ?
    // explicit FaceFieldAll(const MIndex &p)

    FaceFieldAll(const FaceFieldAll &ffc, const MemoryOwner owner)
        : BaseType(ffc, owner)
    {
#ifndef NDEBUG
        assert(this->components_.size() == IndexRangeType::Dim);
        for (auto c : components_) {
            assert(c != nullptr);
        }
#endif /* NDEBUG */
    }

    FaceFieldAll(const std::vector<IndexRangeType> &range_list,
                 const std::vector<DataType *> &ptr_list,
                 const std::vector<size_t> &bytes_list,
                 const std::vector<FieldState *> &state_list)
        : BaseType(range_list, ptr_list, bytes_list, state_list)
    {
#ifndef NDEBUG
        assert(this->components_.size() == IndexRangeType::Dim);
        for (auto c : components_) {
            assert(c != nullptr);
        }
#endif /* NDEBUG */
    }

    FaceFieldAll() = default;
    FaceFieldAll(const FaceFieldAll &c) = default;
    FaceFieldAll(FaceFieldAll &&c) noexcept = default;
    FaceFieldAll &operator=(const FaceFieldAll &c) = default;
    FaceFieldAll &operator=(FaceFieldAll &&c) = default;
    ~FaceFieldAll() = default;

private:
    using BaseType::components_;
};

template <typename TField, size_t RANK>
class TensorField : public FieldContainer<TField>
{
public:
    using BaseType = FieldContainer<TField>;
    using FieldType = typename BaseType::FieldType;
    using FieldState = typename FieldType::FieldStateType;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

private:
    template <size_t B, size_t E>
    struct Power {
        static constexpr size_t value = B * Power<B, E - 1>::value;
    };
    template <size_t B>
    struct Power<B, 0> {
        static constexpr size_t value = 1;
    };
    using MemoryOwner = typename FieldType::BaseType::MemoryOwner;

public:
    static constexpr size_t Rank = RANK;
    static constexpr size_t NComponents =
        Power<IndexRangeType::Dim, RANK>::value;

    explicit TensorField(const IndexRangeType &r)
        : BaseType(NComponents, r, Rank)
    {
    }

    TensorField(const TensorField &ffc, const MemoryOwner owner)
        : BaseType(ffc, owner)
    {
        assert(this->components_.size() == NComponents);
    }

    TensorField(const std::vector<IndexRangeType> &range_list,
                const std::vector<DataType *> &ptr_list,
                const std::vector<size_t> &bytes_list,
                const std::vector<FieldState *> &state_list)
        : BaseType(range_list, ptr_list, bytes_list, state_list)
    {
        assert(this->components_.size() == NComponents);
    }

    TensorField() = default;
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
public:
    using BaseType = TField;
    using FieldType = typename BaseType::FieldType;
    using FieldState = typename FieldType::FieldStateType;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using DataType = typename BaseType::DataType;

private:
    using MemoryOwner = typename FieldType::BaseType::MemoryOwner;

public:
    FieldView(const BaseType &f) : BaseType(f, MemoryOwner::No) {}

    FieldView() = delete;
    FieldView(const FieldView &c) = default;
    FieldView &operator=(const FieldView &c) = default;
    FieldView(FieldView &&c) = delete;
    FieldView &operator=(FieldView &&c) = delete;
    ~FieldView() = default;

    BaseType copy() const { return BaseType(*this, MemoryOwner::Yes); }
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* FIELD_H_7NZ0QFMC */
