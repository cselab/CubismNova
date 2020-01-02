// File       : Field.h
// Created    : Mon Dec 30 2019 05:52:21 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Basic block data field
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef FIELD_H_7NZ0QFMC
#define FIELD_H_7NZ0QFMC

#include "Alloc/AlignedBlockAllocator.h"
#include "Block/Data.h"
#include "Block/FieldOperator.h"
#include "Common.h"

#include <array>
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

/// @brief Meta data (state) of a block field.  Data carried in a field
///        state is supposed to be static.
struct FieldState {
    size_t rank;      // field rank -- rank=0 (scalar), rank=1 (vector), ...
    size_t comp;      // field component in rank dimension
};

// TODO: [fabianw@mavt.ethz.ch; 2020-01-01]
// class FieldUnit // for physical units of carried data

/// @brief Block field base class
///
/// @tparam TBlockData Block data type
/// @tparam TState Field state type
template <typename TBlockData, typename TState = FieldState>
class Field : public TBlockData
{
protected:
    using TBlockData::block_;
    using TBlockData::range_;

    /// @brief Generic iterator for block data
    ///
    /// @tparam T Iterator type
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

    /// @brief Main field constructor
    ///
    /// @param r Index range that spans the data
    /// @param fs Field state (optional)
    explicit Field(const IndexRangeType &r,
                   const FieldStateType &fs = FieldStateType())
        : BaseType(r), state_(new FieldStateType())
    {
        *state_ = fs;
    }

    /// @brief Copy constructor for deep and shallow copies
    ///
    /// @param f Field to be copied
    /// @param o Memory ownership (Data::MemoryOwner::Yes = copy deep)
    Field(const Field &f, const typename BaseType::MemoryOwner o)
        : BaseType(f, o)
    {
        if (this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(f);
    }

    /// @brief Low-level view constructor
    ///
    /// @param r Index range that is spanned by ptr
    /// @param ptr Block data pointer
    /// @param bytes Number of bytes in block data
    /// @param sptr Field state pointer
    Field(const IndexRangeType &r,
          DataType *ptr,
          const size_t bytes,
          FieldStateType *sptr)
        : BaseType(r, ptr, bytes), state_(sptr)
    {
    }

    /// @brief Standard copy constructor
    ///
    /// @param c Field to copy from
    Field(const Field &c) : BaseType(c)
    {
        if (this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(c);
    }

    /// @brief Standard move constructor
    ///
    /// @param c Field to move from
    Field(Field &&c) noexcept
        : BaseType(std::move(c)), state_(std::move(c.state_))
    {
        c.state_ = nullptr;
    }

    /// @brief Virtual destructor
    ~Field() override { disposeState_(); }

    /// @brief Standard copy assignment operator
    ///
    /// @param c Field to copy from
    ///
    /// @return This field with contents copied from c (data and state)
    Field &operator=(const Field &c)
    {
        assert(range_.size() == c.range_.size());
        if (this != &c) {
            copyState_(c);
            BaseType::operator=(c);
        }
        return *this;
    };

    /// @brief Standard move assignment operator
    ///
    /// @param c Field to move from
    ///
    /// @return This field with contents moved from c (data and state)
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

    /// @brief Iterator for underlying data
    using iterator = IteratorBase<DataType>;
    /// @brief Const iterator for underlying data
    using const_iterator = IteratorBase<const DataType>;
    /// @brief Begin of data
    ///
    /// @return Iterator
    iterator begin() noexcept { return iterator(block_); }
    const_iterator begin() const noexcept { return const_iterator(block_); }
    /// @brief End of data
    ///
    /// @return Iterator
    iterator end() noexcept { return iterator(block_ + range_.size()); }
    const_iterator end() const noexcept
    {
        return const_iterator(block_ + range_.size());
    }
    /// @brief Begin of data
    ///
    /// @return Const iterator
    const_iterator cbegin() const noexcept { return const_iterator(block_); }
    /// @brief End of data
    ///
    /// @return Const iterator
    const_iterator cend() const noexcept
    {
        return const_iterator(block_ + range_.size());
    }

    /// @brief Test if field belongs to scalar class
    ///
    /// @return Boolean (true if scalar, false if higher rank tensor class)
    bool isScalar() const { return (0 == state_->rank); }
    /// @brief Rank of field
    ///
    /// @return Rank of tensor field (0 = scalar, 1 = vector, ...)
    size_t getRank() const { return state_->rank; }
    /// @brief Component ID if the field belongs to a set of fields (e.g. higher
    /// rank tensor or field belongs to a field container)
    ///
    /// @return Component ID (integral type >= 0)
    size_t getComp() const { return state_->comp; }
    /// @brief Access to field state
    ///
    /// @return Non-const reference to state
    FieldStateType &getState() { return *state_; }
    /// @brief Access to field state
    ///
    /// @return Const reference to state
    const FieldStateType &getState() const { return *state_; }

    Field operator-() const
    {
        Field f(*this);
        // TODO: [fabianw@mavt.ethz.ch; 2020-01-02] bit-flip might be better
        for (size_t i = 0; i < range_.size(); ++i) {
            f[i] = -f[i];
        }
        return f;
    }

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-02] Field FMA

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

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-02] reciprocal() should not perform
    // the multiplication in fieldRcp.  Factor into specialization
    void reciprocal(const DataType c = 1)
    {
        // reciprocal multiplied by c
        fieldRcp(block_, c, block_, range_.size());
        return;
    }

private:
    FieldStateType *state_; // Field state values/meta data

    /// @brief Deallocate state data
    void disposeState_()
    {
        if (this->isMemoryOwner()) {
            if (state_) {
                delete state_;
            }
        }
        state_ = nullptr;
    }

    /// @brief Copy state data from field c (deep or shallow depending on
    ///  ownership)
    ///
    /// @param c Field to copy state from
    void copyState_(const Field &c)
    {
        if (this->isMemoryOwner()) {
            *state_ = *c.state_;
        } else {
            state_ = c.state_;
        }
    }
};

/// @brief Basic cell-centered data field
///
/// @tparam T Data type (must be POD)
/// @tparam State State type
/// @tparam Alloc Allocator type
template <typename T,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator,
          size_t Dimension = CUBISM_DIMENSION>
using CellField = Field<Data<T, DataMapping::Cell, Dimension, Alloc<T>>, State>;

/// @brief Basic node-centered data field
///
/// @tparam T Data type (must be POD)
/// @tparam State State type
/// @tparam Alloc Allocator type
template <typename T,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator,
          size_t Dimension = CUBISM_DIMENSION>
using NodeField = Field<Data<T, DataMapping::Node, Dimension, Alloc<T>>, State>;

/// @brief Basic face-centered data field.  Faces are stored individually for
/// the dimensionality specified by CUBISM_DIMENSION at compile time.  See
/// the FaceFieldAll type for a container of size CUBISM_DIMENSION.
///
/// @tparam T Data type (must be POD)
/// @tparam State State type
/// @tparam Alloc Allocator type
template <typename T,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator,
          size_t Dimension = CUBISM_DIMENSION>
using FaceField = Field<Data<T, DataMapping::Face, Dimension, Alloc<T>>, State>;

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

/// @brief Field container type (stores pointers to fields).  An incomplete
/// container contains nullptr for some of its components.
///
/// @tparam TField Type of field
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
    /// @brief Default constructor
    ///
    /// @param n Number of components in container
    explicit FieldContainer(const size_t n = 0) : components_(n, nullptr) {}

    /// @brief Construct from a list of pointers
    ///
    /// @param vf Vector of pointers to underlying type (may be nullptr)
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

    /// @brief Copy constructor for deep and shallow copies
    ///
    /// @param fc Field container to be copied
    /// @param o Memory ownership (Data::MemoryOwner::Yes = copy deep)
    FieldContainer(const FieldContainer &fc, const MemoryOwner o)
        : components_(fc.size(), nullptr)
    {
        for (size_t i = 0; i < fc.size(); ++i) {
            const BaseType *comp = fc.components_[i];
            if (comp) {
                components_[i] = new BaseType(*comp, o);
            }
        }
    }

    /// @brief Standard constructor for field container (allocates new fields)
    ///
    /// @param n Number of components
    /// @param r Index range spanned by field data
    /// @param rank Rank of field
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

    /// @brief Low-level constructor for views (never allocates new fields)
    ///
    /// @param range_list Vector of index ranges for each component
    /// @param ptr_list Vector of block data pointer corresponding to range_list
    /// @param bytes_list Number of bytes pointed to by pointer in ptr_list
    /// @param state_list Vector of field state pointers for each component
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

    /// @brief Standard copy constructor
    ///
    /// @param c Field container to copy from
    FieldContainer(const FieldContainer &c) : components_(c.size(), nullptr)
    {
        for (size_t i = 0; i < c.size(); ++i) {
            const BaseType *comp = c.components_[i];
            if (comp) {
                components_[i] = new BaseType(*comp);
            }
        }
    }

    /// @brief Standard move constructor
    ///
    /// @param c Field container to move from
    FieldContainer(FieldContainer &&c) noexcept
        : components_(std::move(c.components_))
    {
        c.components_.clear();
    }

    /// @brief Virtual destructor
    virtual ~FieldContainer() { dispose_(); }

    /// @brief Standard copy assignment operator
    ///
    /// @param rhs Field container to assign from
    ///
    /// @return This field container with copied contents from rhs
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

    /// @brief Standard move assignment operator
    ///
    /// @param rhs Field container to assign from
    ///
    /// @return This field container with moved contents from rhs
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
    /// @brief Iterator for field pointers in ContainerType
    using iterator = typename ContainerType::iterator;
    /// @brief Const iterator for field pointers in ContainerType
    using const_iterator = typename ContainerType::const_iterator;
    /// @brief Reverse iterator for field pointers in ContainerType
    using reverse_iterator = typename ContainerType::reverse_iterator;
    /// @brief Const reverse iterator for field pointers in ContainerType
    using const_reverse_iterator =
        typename ContainerType::const_reverse_iterator;
    /// @brief Begin of fields
    ///
    /// @return Iterator
    iterator begin() noexcept { return components_.begin(); }
    const_iterator begin() const noexcept
    {
        return const_iterator(components_.begin());
    }
    /// @brief End of fields
    ///
    /// @return Iterator
    iterator end() noexcept { return components_.end(); }
    const_iterator end() const noexcept
    {
        return const_iterator(components_.end());
    }
    /// @brief Begin of reversed fields
    ///
    /// @return Reverse iterator
    reverse_iterator rbegin() noexcept { return components_.rbegin(); }
    const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(components_.rbegin());
    }
    /// @brief End of reversed fields
    ///
    /// @return Reverse iterator
    reverse_iterator rend() noexcept { return components_.rend(); }
    const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(components_.rend());
    }
    /// @brief Begin of fields
    ///
    /// @return Const iterator
    const_iterator cbegin() const noexcept { return components_.cbegin(); }
    /// @brief End of fields
    ///
    /// @return Const iterator
    const_iterator cend() const noexcept { return components_.cend(); }
    /// @brief Begin of reversed fields
    ///
    /// @return Const reverse iterator
    const_reverse_iterator crbegin() const noexcept
    {
        return components_.crbegin();
    }
    /// @brief End of reversed fields
    ///
    /// @return Const reverse iterator
    const_reverse_iterator crend() const noexcept
    {
        return components_.crend();
    }

    /// @brief Number of fields contained
    ///
    /// @return Size of container
    size_t size() const { return components_.size(); }

    /// @brief Forced deep copy of underlying fields (to be used with views).
    /// This copies FieldType::BaseType data only, not the state.
    ///
    /// @param rhs Field container to copy from
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

    /// @brief Access to fields.  This method throws a std::runtime_error if the
    /// component is a nullptr.
    ///
    /// @param i Component ID of field to be accessed
    ///
    /// @return Reference to component i
    BaseType &operator[](const size_t i) { return getComp_(i); }
    /// @brief Access to fields.  This method throws a std::runtime_error if the
    /// component is a nullptr.
    ///
    /// @param i Component ID of field to be accessed
    ///
    /// @return Const reference to component i
    const BaseType &operator[](const size_t i) const { return getComp_(i); }

    /// @brief Generic access to fields.  This method throws a
    /// std::runtime_error if the component is a nullptr.
    ///
    /// @tparam T Generic index type that defines casting to size_t
    /// @param t Component ID of field to be accessed
    ///
    /// @return Reference to component t
    template <typename T>
    BaseType &operator[](const T &t)
    {
        return getComp_(static_cast<size_t>(t));
    }
    /// @brief Generic access to fields.  This method throws a
    /// std::runtime_error if the component is a nullptr.
    ///
    /// @tparam T Generic index type that defines casting to size_t
    /// @param t Component ID of field to be accessed
    ///
    /// @return Const reference to component t
    template <typename T>
    const BaseType &operator[](const T &t) const
    {
        return getComp_(static_cast<size_t>(t));
    }

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

    void reciprocal(const DataType c = 1)
    {
        for (size_t i = 0; i < components_.size(); ++i) {
            BaseType *comp = components_[i];
            if (comp) {
                comp->reciprocal(c);
            }
        }
        return;
    }

protected:
    ContainerType components_; // the stars of the show

private:
    /// @brief Return a component
    ///
    /// @param i Component ID
    ///
    /// @return Reference to component
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

    /// @brief Return a component
    ///
    /// @param i Component ID
    ///
    /// @return Const reference to component
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

    /// @brief Deallocation of underlying components
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

/// @brief Container class for all faces in a CUBISM_DIMENSION-ional problem.
/// The underlying face fields are based on the FaceField template.  For
/// CUBISM_DIMENSION \in {1,2,3}, the face field for faces with normal in the X
/// direction can be obtained with ff[0] or ff[Cubism::Dir::X] for example,
/// where ff is of type FaceFieldAll.
///
/// @tparam T Data type of FaceField
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
    /// @brief Main constructor to generate a face field given the cell_domain
    ///
    /// @param cell_domain Index range spanned by the cell domain
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

    /// @brief Copy constructor for deep and shallow copies
    ///
    /// @param ffc Face field container to be copied
    /// @param o Memory ownership (Data::MemoryOwner::Yes = copy deep)
    FaceFieldAll(const FaceFieldAll &ffc, const MemoryOwner o)
        : BaseType(ffc, o)
    {
#ifndef NDEBUG
        assert(this->components_.size() == IndexRangeType::Dim);
        for (auto c : components_) {
            assert(c != nullptr);
        }
#endif /* NDEBUG */
    }

    /// @brief Low-level constructor for views (never allocates new fields)
    ///
    /// @param range_list Vector of index ranges for each component
    /// @param ptr_list Vector of block data pointer corresponding to range_list
    /// @param bytes_list Number of bytes pointed to by pointer in ptr_list
    /// @param state_list Vector of field state pointers for each component
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

    /// @brief Default constructor generates an empty container
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

    /// @brief Main constructor to generate a tensor field
    ///
    /// @param r Index range spanned by the field data
    explicit TensorField(const IndexRangeType &r)
        : BaseType(NComponents, r, Rank)
    {
    }

    /// @brief Copy constructor for deep and shallow copies
    ///
    /// @param tfc Tensor field container to be copied
    /// @param o Memory ownership (Data::MemoryOwner::Yes = copy deep)
    TensorField(const TensorField &tfc, const MemoryOwner o) : BaseType(tfc, o)
    {
        assert(this->components_.size() == NComponents);
    }

    /// @brief Low-level constructor for views (never allocates new fields)
    ///
    /// @param range_list Vector of index ranges for each component
    /// @param ptr_list Vector of block data pointer corresponding to range_list
    /// @param bytes_list Number of bytes pointed to by pointer in ptr_list
    /// @param state_list Vector of field state pointers for each component
    TensorField(const std::vector<IndexRangeType> &range_list,
                const std::vector<DataType *> &ptr_list,
                const std::vector<size_t> &bytes_list,
                const std::vector<FieldState *> &state_list)
        : BaseType(range_list, ptr_list, bytes_list, state_list)
    {
        assert(this->components_.size() == NComponents);
    }

    /// @brief Default constructor generates an empty container
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

/// @brief Convenience type for vector fields
///
/// @tparam TField Field type
template <typename TField>
using VectorField = TensorField<TField, 1>;

/// @brief Field view type (never owns memory)
///
/// @tparam TField Field type
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
    /// @brief Main constructor to generate a field view
    ///
    /// @param f Field for which to generate the view
    FieldView(const BaseType &f) : BaseType(f, MemoryOwner::No) {}

    FieldView() = delete;
    FieldView(const FieldView &c) = default;
    FieldView &operator=(const FieldView &c) = default;
    FieldView(FieldView &&c) = delete;
    FieldView &operator=(FieldView &&c) = delete;
    ~FieldView() = default;

    /// @brief Force a deep copy of the viewed field
    ///
    /// @return Copy (new memory) of field view
    BaseType copy() const { return BaseType(*this, MemoryOwner::Yes); }
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* FIELD_H_7NZ0QFMC */
