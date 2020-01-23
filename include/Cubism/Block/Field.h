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

#ifndef NDEBUG
#include <cstdio>
#endif /* NDEBUG */

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Block)

/**
 * @brief Default meta data (state) of a block field
 *
 * @rst
 * Minimal (default) state of a field is empty.  Custom field state types my add
 * additional state to describe meta data of a field.  A field state type must
 * define copy semantics.  Components in ``TensorField`` types and
 * ``FaceContainer`` types share one instance of a state.
 * @endrst
 */
struct FieldState {
};

// TODO: [fabianw@mavt.ethz.ch; 2020-01-01]
// class FieldUnit // for physical units of carried data

/**
 * @brief Block field base class
 * @tparam TBlockData Block data type
 * @tparam TState Field state type
 *
 * @rst
 * Generic block field type used by :ref:`grid` classes to compose a certain
 * topology of block fields.
 * @endrst
 */
template <typename TBlockData, typename TState = FieldState>
class Field : public TBlockData
{
protected:
    using TBlockData::block_;
    using TBlockData::range_;

    /**
     * @brief Generic iterator for block data
     * @tparam T Iterator type
     */
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
    using BlockDataType = TBlockData;
    using typename BaseType::DataType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using FieldStateType = TState;

    static constexpr size_t Rank = 0;
    static constexpr size_t NComponents = 1;
    static constexpr Cubism::EntityType EntityType = BlockDataType::EntityType;

    Field() = delete;

    /**
     * @brief Main field constructor
     * @param r Index range that spans the data
     * @param fs Field state
     */
    explicit Field(const IndexRangeType &r,
                   const FieldStateType &fs = FieldStateType())
        : BaseType(r), state_(new FieldStateType())
    {
        *state_ = fs;
    }

    /**
     * @brief Copy constructor for deep and shallow copies
     * @param f Field to be copied
     * @param o Memory ownership (Data::MemoryOwner::Yes = copy deep)
     */
    Field(const Field &f, const typename BaseType::MemoryOwner o)
        : BaseType(f, o)
    {
        if (this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(f);
    }

    /**
     * @brief Low-level constructor for external memory
     * @param r Index range that is spanned by ptr
     * @param ptr Block data pointer
     * @param bytes Number of bytes in block data
     * @param sptr Field state pointer
     */
    Field(const IndexRangeType &r,
          DataType *ptr,
          const size_t bytes,
          FieldStateType *sptr)
        : BaseType(r, ptr, bytes), state_(sptr)
    {
    }

    /**
     * @brief Low-level constructor for external memory
     * @param range_list Vector of index range
     * @param ptr_list Vector of block data pointer
     * @param bytes_list Number of bytes pointed to by pointer in ptr_list
     * @param state_list Vector of field state pointer
     */
    Field(const std::vector<IndexRangeType> &range_list,
          const std::vector<DataType *> &ptr_list,
          const std::vector<size_t> &bytes_list,
          const std::vector<FieldStateType *> &state_list)
        : BaseType(range_list[0], ptr_list[0], bytes_list[0]),
          state_(state_list[0])
    {
        assert(range_list.size() == 1);
    }

    /**
     * @brief Standard copy constructor
     * @param c Field to copy from
     */
    Field(const Field &c) : BaseType(c)
    {
        if (this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(c);
    }

    /**
     * @brief Standard move constructor
     * @param c Field to move from
     */
    Field(Field &&c) noexcept
        : BaseType(std::move(c)), state_(std::move(c.state_))
    {
        c.state_ = nullptr;
    }

    /**
     * @brief Virtual destructor
     */
    ~Field() override
    {
        if (this->external_memory_) {
            // state is not deallocated in the case of externally managed memory
            state_ = nullptr;
        }
        disposeState_();
    }

    /**
     * @brief Standard copy assignment operator
     * @param c Field to copy from
     * @return This field with contents copied from c (data and state)
     */
    Field &operator=(const Field &c)
    {
        assert(range_.size() == c.range_.size());
        if (this != &c) {
            copyState_(c);
            BaseType::operator=(c);
        }
        return *this;
    };

    /**
     * @brief Standard move assignment operator
     * @param c Field to move from
     * @return This field with contents moved from c (data and state)
     */
    Field &operator=(Field &&c)
    {
        assert(range_.size() == c.range_.size());
        if (this != &c) {
            BaseType::operator=(std::move(c));
            if (this->external_memory_) {
                // state is not deallocated in the case of externally managed
                // memory
                state_ = nullptr;
            }
            disposeState_();
            state_ = std::move(c.state_);
            c.state_ = nullptr;
        }
        return *this;
    }

    /**
     * @brief Iterator for underlying data
     */
    using iterator = IteratorBase<DataType>;
    /**
     * @brief Const iterator for underlying data
     */
    using const_iterator = IteratorBase<const DataType>;
    /**
     * @brief Begin of data
     * @return Iterator
     */
    iterator begin() noexcept { return iterator(block_); }
    const_iterator begin() const noexcept { return const_iterator(block_); }
    /**
     * @brief End of data
     * @return Iterator
     */
    iterator end() noexcept { return iterator(block_ + range_.size()); }
    const_iterator end() const noexcept
    {
        return const_iterator(block_ + range_.size());
    }
    /**
     * @brief Begin of data
     * @return Const iterator
     */
    const_iterator cbegin() const noexcept { return const_iterator(block_); }
    /**
     * @brief End of data
     * @return Const iterator
     */
    const_iterator cend() const noexcept
    {
        return const_iterator(block_ + range_.size());
    }

    /**
     * @brief Number of data elements carried by field
     * @return Size of data index range
     */
    size_t size() const { return range_.size(); }

    /**
     * @brief Test if field belongs to scalar class
     * @return Boolean (true if scalar, false if higher rank tensor class)
     */
    bool isScalar() const { return (0 == Rank); }
    /**
     * @brief Access to field state
     * @return Non-const reference to state
     */
    FieldStateType &getState() { return *state_; }
    /**
     * @brief Access to field state
     * @return ``const`` reference to state
     */
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

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-17] document these members
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

    /**
     * @brief Deallocate state data
     */
    void disposeState_()
    {
        if (this->isMemoryOwner()) {
            if (state_) {
                delete state_;
            }
        }
        state_ = nullptr;
    }

    /**
     * @brief Copy state data from field c (deep or shallow depending on
     *  ownership)
     * @param c Field to copy state from
     */
    void copyState_(const Field &c)
    {
        if (this->isMemoryOwner()) {
            *state_ = *c.state_;
        } else {
            state_ = c.state_;
        }
    }
};

template <typename TBlockData, typename TState>
constexpr size_t Field<TBlockData, TState>::Rank;

template <typename TBlockData, typename TState>
constexpr size_t Field<TBlockData, TState>::NComponents;

template <typename TBlockData, typename TState>
constexpr Cubism::EntityType Field<TBlockData, TState>::EntityType;

/**
 * @brief Basic cell-centered data field
 * @tparam T Data type (must be POD)
 * @tparam State State type
 * @tparam Dimension
 * @tparam Alloc Allocator type
 */
template <typename T,
          typename State = FieldState,
          size_t Dimension = CUBISM_DIMENSION,
          template <typename> class Alloc = AlignedBlockAllocator>
using CellField = Field<Data<T, EntityType::Cell, Dimension, Alloc<T>>, State>;

/**
 * @brief Basic node-centered data field
 * @tparam T Data type (must be POD)
 * @tparam State State type
 * @tparam Dimension
 * @tparam Alloc Allocator type
 */
template <typename T,
          typename State = FieldState,
          size_t Dimension = CUBISM_DIMENSION,
          template <typename> class Alloc = AlignedBlockAllocator>
using NodeField = Field<Data<T, EntityType::Node, Dimension, Alloc<T>>, State>;

/**
 * @brief Basic face-centered data field.
 * @tparam T Data type (must be POD)
 * @tparam State State type
 * @tparam Dimension
 * @tparam Alloc Allocator type
 *
 * @rst
 * Faces are stored individually for the dimensionality specified by
 * ``CUBISM_DIMENSION`` at compile time.  See the ``FaceContainer`` type for a
 * container of size ``CUBISM_DIMENSION``.
 * @endrst
 */
template <typename T,
          typename State = FieldState,
          size_t Dimension = CUBISM_DIMENSION,
          template <typename> class Alloc = AlignedBlockAllocator>
using FaceField = Field<Data<T, EntityType::Face, Dimension, Alloc<T>>, State>;

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

/**
 * @brief Field container type.
 * @tparam TField Type of field
 *
 * @rst
 * This is an actively managed field container.  Unlike ``std::vector``, the
 * destructors of the container elements are called upon destruction of the
 * container.  An incomplete container contains ``nullptr`` for some of its
 * components.
 * @endrst
 */
template <typename TField>
class FieldContainer
{
public:
    using BaseType = TField;
    using DataType = typename BaseType::DataType;
    using FieldType = typename BaseType::FieldType;
    using BlockDataType = typename FieldType::BlockDataType;
    using IndexRangeType = typename BaseType::IndexRangeType;
    using MultiIndex = typename BaseType::MultiIndex;
    using MemoryOwner = typename TField::BaseType::MemoryOwner;

private:
    using ContainerType = std::vector<BaseType *>;
    using FieldStateType = typename FieldType::FieldStateType;

public:
    /**
     * @brief Construct from a list of pointers
     * @param vf Vector of pointers to underlying type (may be ``nullptr``)
     */
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

    /**
     * @brief Copy constructor for deep and shallow copies
     * @param fc Field container to be copied
     * @param o Memory ownership (``Data::MemoryOwner::Yes`` = copy deep)
     */
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

    /**
     * @brief Standard constructor for field container (allocates new fields)
     * @param n Number of components
     * @param r Index range spanned by field data
     * @param rank Rank of field
     */
    FieldContainer(const size_t n,
                   const IndexRangeType &r,
                   const size_t rank = 0)
        : components_(n, nullptr)
    {
        FieldStateType fs;
        fs.rank = rank;
        for (size_t i = 0; i < n; ++i) {
            fs.comp = i;
            components_[i] = new FieldType(r, fs);
        }
    }

    /**
     * @brief Low-level constructor for external memory
     * @param r Index range that is spanned by ``ptr``
     * @param ptr Block data pointer
     * @param bytes Number of bytes in block data
     * @param sptr Field state pointer
     * @param size Size of container
     */
    FieldContainer(const IndexRangeType &r,
                   DataType *ptr,
                   const size_t bytes,
                   FieldStateType *sptr,
                   const size_t size)
        : components_(size, nullptr)
    {
        for (size_t i = 0; i < components_.size(); ++i) {
            assert(ptr != nullptr);
            components_[i] = new FieldType(r, ptr, bytes, sptr);
        }
    }

    /**
     * @brief Low-level constructor for external memory
     * @param range_list Vector of index ranges for each component
     * @param ptr_list Vector of block data pointer corresponding to
     * ``range_list``
     * @param bytes_list Number of bytes pointed to by pointer in ``ptr_list``
     * @param state_list Vector of field state pointers for each component
     */
    FieldContainer(const std::vector<IndexRangeType> &range_list,
                   const std::vector<DataType *> &ptr_list,
                   const std::vector<size_t> &bytes_list,
                   const std::vector<FieldStateType *> &state_list)
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

    /**
     * @brief Standard copy constructor
     * @param c Field container to copy from
     */
    FieldContainer(const FieldContainer &c) : components_(c.size(), nullptr)
    {
        for (size_t i = 0; i < c.size(); ++i) {
            const BaseType *comp = c.components_[i];
            if (comp) {
                components_[i] = new BaseType(*comp);
            }
        }
    }

    /**
     * @brief Standard move constructor
     * @param c Field container to move from
     */
    FieldContainer(FieldContainer &&c) noexcept
        : components_(std::move(c.components_))
    {
        c.components_.clear();
    }

    /**
     * @brief Default constructor
     */
    FieldContainer() = default;

    /**
     * @brief Virtual destructor
     */
    virtual ~FieldContainer() { dispose_(); }

    /**
     * @brief Standard copy assignment operator
     * @param rhs Field container to assign from
     * @return This field container with copied contents from ``rhs``
     */
    FieldContainer &operator=(const FieldContainer &rhs)
    {
        if (this != &rhs) {
            if (components_.size() != rhs.components_.size()) {
#ifndef NDEBUG
                if (components_.size() > 0) {
                    std::fprintf(
                        stderr,
                        "DEBUG: You are using FieldContainer::operator=() with "
                        "two containers of unequal size and the destination "
                        "size is not zero. Is this your intention?\n");
                }
#endif /* NDEBUG */
                dispose_();
                components_.resize(rhs.size());
                for (size_t i = 0; i < rhs.size(); ++i) {
                    components_[i] = nullptr;
                    BaseType *comp = rhs.components_[i];
                    if (comp) {
                        components_[i] = new BaseType(*comp);
                    }
                }
            } else {
                for (size_t i = 0; i < rhs.size(); ++i) {
                    assert(components_[i] != nullptr &&
                           rhs.components_[i] != nullptr);
                    *components_[i] = *rhs.components_[i];
                }
            }
        }
        return *this;
    }

    /**
     * @brief Standard move assignment operator
     * @param rhs Field container to assign from
     * @return This field container with moved contents from ``rhs``
     */
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
    /**
     * @brief Iterator for field pointers in ContainerType
     */
    using iterator = typename ContainerType::iterator;
    /**
     * @brief ``const`` iterator for field pointers in ContainerType
     */
    using const_iterator = typename ContainerType::const_iterator;
    /**
     * @brief Reverse iterator for field pointers in ContainerType
     */
    using reverse_iterator = typename ContainerType::reverse_iterator;
    /**
     * @brief ``const`` reverse iterator for field pointers in ContainerType
     */
    using const_reverse_iterator =
        typename ContainerType::const_reverse_iterator;
    /**
     * @brief Begin of fields
     * @return Iterator
     */
    iterator begin() noexcept { return components_.begin(); }
    const_iterator begin() const noexcept
    {
        return const_iterator(components_.begin());
    }
    /**
     * @brief End of fields
     * @return Iterator
     */
    iterator end() noexcept { return components_.end(); }
    const_iterator end() const noexcept
    {
        return const_iterator(components_.end());
    }
    /**
     * @brief Begin of reversed fields
     * @return Reverse iterator
     */
    reverse_iterator rbegin() noexcept { return components_.rbegin(); }
    const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(components_.rbegin());
    }
    /**
     * @brief End of reversed fields
     * @return Reverse iterator
     */
    reverse_iterator rend() noexcept { return components_.rend(); }
    const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(components_.rend());
    }
    /**
     * @brief Begin of fields
     * @return ``const`` iterator
     */
    const_iterator cbegin() const noexcept { return components_.cbegin(); }
    /**
     * @brief End of fields
     * @return ``const`` iterator
     */
    const_iterator cend() const noexcept { return components_.cend(); }
    /**
     * @brief Begin of reversed fields
     * @return ``const`` reverse iterator
     */
    const_reverse_iterator crbegin() const noexcept
    {
        return components_.crbegin();
    }
    /**
     * @brief End of reversed fields
     * @return ``const`` reverse iterator
     */
    const_reverse_iterator crend() const noexcept
    {
        return components_.crend();
    }

    /**
     * @brief Number of fields contained
     * @return Size of container
     */
    size_t size() const { return components_.size(); }

    /**
     * @brief Forced deep copy of underlying fields
     * @param rhs Field container to copy from
     *
     * @rst
     * This copies ``FieldType::BaseType`` data only.
     *
     * Example: ``fv`` is a field view and ``f`` is another field (either view
     * or memory owner)
     *
     * .. code-block:: cpp
     *
     *    fv = f; // shallow copy of f to fv [updates pointers in fv only]
     *    fv.copyData(f); // deep copy of data in f to fv [expensive operation]
     * @endrst
     */
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

    /**
     * @brief Append new component at the end of container
     * @param p Pointer to new component
     *
     * @rst
     * The destructor of ``FieldContainer`` calls the destructor of ``p`` when
     * the object is destroyed.
     * @endrst
     */
    void pushBack(BaseType *p) { components_.push_back(p); }

    /**
     * @brief Destroy container components
     *
     * @rst
     * Unlike to ``std::vector``, this calls the destructor of each container
     * component effectively destroying the object.  The size of the container
     * after this operation is zero.
     * @endrst
     */
    void clear()
    {
        dispose_();
        this->components_.clear();
    }

    /**
     * @brief Access to fields
     * @param i Component ID of field to be accessed
     * @return Reference to component ``i``
     *
     * @rst
     * This method throws a ``std::runtime_error`` if the component is a
     * ``nullptr``.
     * @endrst
     */
    BaseType &operator[](const size_t i) { return getComp_(i); }
    /**
     * @brief Access to fields
     * @param i Component ID of field to be accessed
     * @return ``const`` reference to component ``i``
     *
     * @rst
     * This method throws a ``std::runtime_error`` if the component is a
     * ``nullptr``.
     * @endrst
     */
    const BaseType &operator[](const size_t i) const { return getComp_(i); }

    /**
     * @brief Generic access to fields.
     * @tparam T Generic index type that defines casting to ``size_t``
     * @param t Component ID of field to be accessed
     * @return Reference to component ``t``
     *
     * @rst
     * This method throws a ``std::runtime_error`` if the component is a
     * ``nullptr``.
     * @endrst
     */
    template <typename T>
    BaseType &operator[](const T &t)
    {
        return getComp_(static_cast<size_t>(t));
    }
    /**
     * @brief Generic access to fields.
     * @tparam T Generic index type that defines casting to ``size_t``
     * @param t Component ID of field to be accessed
     * @return ``const`` reference to component ``t``
     *
     * @rst
     * This method throws a ``std::runtime_error`` if the component is a
     * ``nullptr``.
     * @endrst
     */
    template <typename T>
    const BaseType &operator[](const T &t) const
    {
        return getComp_(static_cast<size_t>(t));
    }

    // TODO: [fabianw@mavt.ethz.ch; 2020-01-17] document these methods
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
    /**
     * @brief Return a component
     * @param i Component ID
     * @return Reference to component
     */
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

    /**
     * @brief Return a component
     * @param i Component ID
     * @return ``const`` reference to component
     */
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

    /**
     * @brief Deallocation of underlying components
     */
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

/**
 * @brief Container class for all faces in a ``CUBISM_DIMENSION``-ional problem
 * @tparam T Data type of underlying ``FaceField``
 *
 * @rst
 * The underlying face fields are based on the ``FaceField`` template.  For
 * ``CUBISM_DIMENSION`` in ``{1,2,3}``, the face field for faces with normal in
 * the ``X`` direction can be obtained with ``ff[0]`` or ``ff[Cubism::Dir::X]``
 * for example, where ``ff`` is of type ``FaceContainer``.
 * @endrst
 */
template <typename T,
          typename State = FieldState,
          size_t Dimension = CUBISM_DIMENSION,
          template <typename> class Alloc = AlignedBlockAllocator>
class FaceContainer
    : public FieldContainer<FaceField<T, State, Dimension, Alloc>>
{
public:
    using BaseType = FieldContainer<FaceField<T, State, Dimension, Alloc>>;
    using typename BaseType::BlockDataType;
    using typename BaseType::DataType;
    using typename BaseType::FieldType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using FieldStateType = typename FieldType::FieldStateType;
    using MemoryOwner = typename FieldType::BaseType::MemoryOwner;

    static constexpr size_t Rank = FieldType::Rank;
    static constexpr size_t NComponents = FieldType::NComponents;
    static constexpr Cubism::EntityType EntityType = BlockDataType::EntityType;

    /**
     * @brief Main constructor to generate a face field given the
     * ``cell_domain``
     * @param cell_domain Index range spanned by the cell domain
     */
    explicit FaceContainer(const IndexRangeType &cell_domain)
    {
        const MultiIndex cells = cell_domain.getExtent(); // number of cells
        FieldStateType fs;
        fs.rank = 0;
        for (size_t i = 0; i < IndexRangeType::Dim; ++i) {
            // XXX: [fabianw@mavt.ethz.ch; 2020-01-01] Not most favorable for
            // vectorization but more intuitive.  Will require some
            // transposition in vectorized code.
            const IndexRangeType ri(cells + MultiIndex::getUnitVector(i));
            fs.comp = i;
            components_.push_back(new FieldType(ri, fs));
            assert(components_[i] != nullptr);
        }
    }

    /**
     * @brief Copy constructor for deep and shallow copies
     * @param ffc Face field container to be copied
     * @param o Memory ownership (``Data::MemoryOwner::Yes`` = copy deep)
     */
    FaceContainer(const FaceContainer &ffc, const MemoryOwner o)
        : BaseType(ffc, o)
    {
#ifndef NDEBUG
        assert(this->components_.size() == IndexRangeType::Dim);
        for (auto c : components_) {
            assert(c != nullptr);
        }
#endif /* NDEBUG */
    }

    /**
     * @brief Low-level constructor for external memory
     * @param r Index range that is spanned by ``ptr``
     * @param ptr Block data pointer
     * @param bytes Number of bytes in block data
     * @param sptr Field state pointer
     */
    FaceContainer(const IndexRangeType &r,
                  DataType *ptr,
                  const size_t bytes,
                  FieldStateType *sptr)
        : BaseType(r, ptr, bytes, sptr, IndexRangeType::Dim)
    {
#ifndef NDEBUG
        assert(this->components_.size() == IndexRangeType::Dim);
        for (auto c : components_) {
            assert(c != nullptr);
        }
#endif /* NDEBUG */
    }

    /**
     * @brief Low-level constructor for external memory
     * @param range_list Vector of index ranges for each component
     * @param ptr_list Vector of block data pointer corresponding to
     * ``range_list``
     * @param bytes_list Number of bytes pointed to by pointer in ``ptr_list``
     * @param state_list Vector of field state pointers for each component
     */
    FaceContainer(const std::vector<IndexRangeType> &range_list,
                  const std::vector<DataType *> &ptr_list,
                  const std::vector<size_t> &bytes_list,
                  const std::vector<FieldStateType *> &state_list)
        : BaseType(range_list, ptr_list, bytes_list, state_list)
    {
#ifndef NDEBUG
        assert(this->components_.size() == IndexRangeType::Dim);
        for (auto c : components_) {
            assert(c != nullptr);
        }
#endif /* NDEBUG */
    }

    /**
     * @brief Default constructor generates an empty container
     */
    FaceContainer() = default;
    FaceContainer(const FaceContainer &c) = default;
    FaceContainer(FaceContainer &&c) = default;
    FaceContainer &operator=(const FaceContainer &c) = default;
    FaceContainer &operator=(FaceContainer &&c) = default;
    ~FaceContainer() = default;

    /**
     * @brief Face field access
     * @param i Direction index
     */
    FieldType &operator[](const size_t i)
    {
        assert(i < IndexRangeType::Dim);
        return *components_[i];
    }

    /**
     * @brief Face field access
     * @param i Direction index
     */
    const FieldType &operator[](const size_t i) const
    {
        assert(i < IndexRangeType::Dim);
        return *components_[i];
    }

    /**
     * @brief Face field access
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Direction of face
     */
    template <typename Dir>
    FieldType &operator[](const Dir &t)
    {
        return this->operator[](static_cast<size_t>(t));
    }

    /**
     * @brief Face field access
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Direction of face
     */
    template <typename Dir>
    const FieldType &operator[](const Dir &t) const
    {
        return this->operator[](static_cast<size_t>(t));
    }

private:
    using BaseType::components_;
};

template <typename T,
          typename State,
          size_t Dimension,
          template <typename>
          class Alloc>
constexpr size_t FaceContainer<T, State, Dimension, Alloc>::Rank;

template <typename T,
          typename State,
          size_t Dimension,
          template <typename>
          class Alloc>
constexpr size_t FaceContainer<T, State, Dimension, Alloc>::NComponents;

template <typename T,
          typename State,
          size_t Dimension,
          template <typename>
          class Alloc>
constexpr Cubism::EntityType
    FaceContainer<T, State, Dimension, Alloc>::EntityType;

/**
 * @brief Generic tensor field
 * @tparam TField Field type
 * @tparam RANK Tensor rank
 */
template <typename TField, size_t RANK>
class TensorField : public FieldContainer<TField>
{
public:
    using BaseType = FieldContainer<TField>;
    using TensorComponentType = TField;
    using typename BaseType::BlockDataType;
    using typename BaseType::DataType;
    using typename BaseType::FieldType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using FieldStateType = typename FieldType::FieldStateType;
    using MemoryOwner = typename FieldType::BaseType::MemoryOwner;

private:
    template <size_t B, size_t E>
    struct Power {
        static constexpr size_t value = B * Power<B, E - 1>::value;
    };
    template <size_t B>
    struct Power<B, 0> {
        static constexpr size_t value = 1;
    };
    using BaseType::components_;

public:
    static constexpr size_t Rank = RANK;
    static constexpr size_t NComponents =
        Power<IndexRangeType::Dim, RANK>::value;
    static constexpr Cubism::EntityType EntityType = BlockDataType::EntityType;

    /**
     * @brief Main constructor to generate a tensor field
     * @param r Index range spanned by the field data
     */
    explicit TensorField(const IndexRangeType &r)
        : BaseType(NComponents, r, Rank)
    {
    }

    /**
     * @brief Copy constructor for deep and shallow copies
     * @param tfc Tensor field container to be copied
     * @param o Memory ownership (``Data::MemoryOwner::Yes`` = copy deep)
     */
    TensorField(const TensorField &tfc, const MemoryOwner o) : BaseType(tfc, o)
    {
        assert(this->components_.size() == NComponents);
    }

    /**
     * @brief Low-level constructor for external memory
     * @param range_list Vector of index ranges for each component
     * @param ptr_list Vector of block data pointer corresponding to
     * ``range_list``
     * @param bytes_list Number of bytes pointed to by pointer in ``ptr_list``
     * @param state_list Vector of field state pointers for each component
     */
    TensorField(const std::vector<IndexRangeType> &range_list,
                const std::vector<DataType *> &ptr_list,
                const std::vector<size_t> &bytes_list,
                const std::vector<FieldStateType *> &state_list)
        : BaseType(range_list, ptr_list, bytes_list, state_list)
    {
        assert(this->components_.size() == NComponents);
    }

    /**
     * @brief Default constructor generates an empty container
     */
    TensorField() = default;
    TensorField(const TensorField &c) = default;
    TensorField(TensorField &&c) = default;
    TensorField &operator=(const TensorField &c) = default;
    TensorField &operator=(TensorField &&c) = default;
    ~TensorField() = default;

    /**
     * @brief Component field access
     * @param i Component index
     */
    FieldType &operator[](const size_t i)
    {
        assert(i < NComponents);
        return *components_[i];
    }

    /**
     * @brief Component field access
     * @param i Component index
     */
    const FieldType &operator[](const size_t i) const
    {
        assert(i < NComponents);
        return *components_[i];
    }

    /**
     * @brief Component field access
     * @tparam TComp Special type that defines a cast to ``size_t``
     * @param t Component of tensor
     */
    template <typename TComp>
    FieldType &operator[](const TComp &t)
    {
        return this->operator[](static_cast<size_t>(t));
    }

    /**
     * @brief Component field access
     * @tparam TComp Special type that defines a cast to ``size_t``
     * @param t Component of tensor
     */
    template <typename TComp>
    const FieldType &operator[](const TComp &t) const
    {
        return this->operator[](static_cast<size_t>(t));
    }
};

template <typename TField, size_t RANK>
constexpr size_t TensorField<TField, RANK>::Rank;

template <typename TField, size_t RANK>
constexpr size_t TensorField<TField, RANK>::NComponents;

template <typename TField, size_t RANK>
constexpr Cubism::EntityType TensorField<TField, RANK>::EntityType;

/**
 * @brief Convenience type for vector fields
 * @tparam TField Field type
 */
template <typename TField>
using VectorField = TensorField<TField, 1>;

/**
 * @brief Field view type
 * @tparam TField Field type
 */
template <typename TField>
class FieldView : public TField
{
public:
    using BaseType = TField;
    using typename BaseType::DataType;
    using typename BaseType::FieldType;
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using FieldStateType = typename FieldType::FieldStateType;
    using MemoryOwner = typename FieldType::BaseType::MemoryOwner;

    /**
     * @brief Main constructor to generate a field view
     * @param f Field for which to generate the view
     */
    FieldView(const BaseType &f) : BaseType(f, MemoryOwner::No) {}

    FieldView() = delete;
    FieldView(const FieldView &c) = default;
    FieldView &operator=(const FieldView &c) = default;
    FieldView(FieldView &&c) = delete;
    FieldView &operator=(FieldView &&c) = delete;
    ~FieldView() = default;

    /**
     * @brief Set new internal view
     * @param c Base field be be viewed at
     */
    void setView(const BaseType &c) { BaseType::operator=(c); }

    /**
     * @brief Force a deep copy of the viewed field
     * @return Deep copy (new allocation) of field view
     */
    BaseType copy() const { return BaseType(*this, MemoryOwner::Yes); }
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* FIELD_H_7NZ0QFMC */
