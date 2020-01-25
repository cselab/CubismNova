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
#include <stdexcept>
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

/** @brief Block scalar field base class
 * @tparam T Field data type
 * @tparam Entity Entity type
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator
 *
 * @rst
 * Generic block scalar field type used by :ref:`grid` classes to compose a
 * certain topology of block fields.
 * @endrst
 * */
template <typename T,
          Cubism::EntityType Entity,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
class Field : public Data<T, Entity, DIM, Alloc<T>>
{
public:
    using FieldType = Field; // scalar field
    using BlockDataType = Data<T, Entity, DIM, Alloc<T>>;
    using typename BlockDataType::DataType;
    using typename BlockDataType::IndexRangeType;
    using typename BlockDataType::MultiIndex;
    using FieldStateType = State;

protected:
    template <typename U>
    using NestedVector = std::vector<std::vector<U>>; // simplifies construction

    using BlockDataType::block_;
    using BlockDataType::range_;

    /**
     * @brief Generic iterator for block data
     * @tparam U Iterator base type
     */
    template <typename U>
    class IteratorBase
    {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = U;
        using difference_type = std::ptrdiff_t;
        using pointer = U *;
        using reference = U &;

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
    static constexpr size_t Rank = 0;
    static constexpr size_t NComponents = 1;
    static constexpr Cubism::EntityType EntityType = BlockDataType::EntityType;
    static constexpr Cubism::FieldClass Class = Cubism::FieldClass::Scalar;
    static_assert(IndexRangeType::Dim > 0, "DIM must be greater than zero");

    /** @brief Default constructor */
    Field() = delete;

    /**
     * @brief Main field constructor
     * @param r Index range that spans the data
     * @param fs Field state
     */
    explicit Field(const IndexRangeType &r,
                   const FieldStateType &fs = FieldStateType())
        : BlockDataType(r), is_subfield_(false), state_(new FieldStateType())
    {
        *state_ = fs;
    }

    /**
     * @brief Low-level field constructor (for higher rank tensors)
     * @param r Index range that spans the data
     * @param pfs Field state pointer (externally managed)
     */
    explicit Field(const IndexRangeType &r, FieldStateType *pfs)
        : BlockDataType(r), is_subfield_(true), state_(pfs)
    {
    }

    /**
     * @brief Low-level copy constructor for deep and shallow copies
     * @param f Field to be copied
     * @param o Memory ownership (``BlockDataType::MemoryOwner::Yes``: deep
     * copy)
     */
    Field(const Field &f, const typename BlockDataType::MemoryOwner o)
        : BlockDataType(f, o), is_subfield_(false), state_(nullptr)
    {
        if (!is_subfield_ && this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(f);
    }

    /**
     * @brief Low-level copy constructor for deep and shallow copies (for higher
     * rank tensors)
     * @param f Field to be copied
     * @param o Memory ownership (``BlockDataType::MemoryOwner::Yes``: deep
     * copy)
     * @param pfs External field state
     */
    Field(const Field &f,
          const typename BlockDataType::MemoryOwner o,
          FieldStateType *pfs)
        : BlockDataType(f, o), is_subfield_(true), state_(pfs)
    {
    }

    /**
     * @brief Low-level constructor for external memory management
     * @param range_list Vector of index range
     * @param ptr_list Vector of block data pointer
     * @param bytes_list Number of bytes pointed to by pointer in ptr_list
     * @param state_list Vector of field state pointer
     * @param subfield Subfield component indicator
     *
     * The nested vector data structure is used to simplify the interface
     * between scalar fields, tensor fields and face field containers.  A
     * subfield component applies to higher rank tensors and face field
     * containers which share the field state among each others.
     */
    Field(const NestedVector<IndexRangeType> &range_list,
          const NestedVector<DataType *> &ptr_list,
          const NestedVector<size_t> &bytes_list,
          const NestedVector<FieldStateType *> &state_list,
          const bool subfield = false)
        : BlockDataType(range_list[0][0], ptr_list[0][0], bytes_list[0][0]),
          is_subfield_(subfield), state_(state_list[0][0])
    {
        // outer
        assert(range_list.size() == 1);
        assert(ptr_list.size() == 1);
        assert(bytes_list.size() == 1);
        assert(state_list.size() == 1);

        // inner
        assert(range_list[0].size() == 1);
        assert(ptr_list[0].size() == 1);
        assert(bytes_list[0].size() == 1);
        assert(state_list[0].size() == 1);
    }

    /**
     * @brief Standard copy constructor for a scalar field
     * @param c Field to copy from
     *
     * This constructor is not designed to be used with individual components of
     * a rank > 0 tensor or face containers.  The copy constructors of these
     * data structures should be used instead (unless you know what you do).
     */
    Field(const Field &c)
        : BlockDataType(c), is_subfield_(c.is_subfield_), state_(nullptr)
    {
        if (!is_subfield_ && this->isMemoryOwner()) {
            state_ = new FieldStateType();
        }
        copyState_(c);
    }

    /**
     * @brief Standard move constructor
     * @param c Field to move from
     */
    Field(Field &&c) noexcept
        : BlockDataType(std::move(c)), is_subfield_(c.is_subfield_),
          state_(std::move(c.state_))
    {
        c.state_ = nullptr;
    }

    /**
     * @brief Virtual destructor
     */
    ~Field() override
    {
        if (is_subfield_ || !this->isMemoryOwner() || this->external_memory_) {
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
            if (!is_subfield_) {
                copyState_(c);
            } else if (is_subfield_ && !this->isMemoryOwner()) {
                state_ = c.state_;
            }
            BlockDataType::operator=(c);
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
            BlockDataType::operator=(std::move(c));
            if (is_subfield_ || !this->isMemoryOwner()) {
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
     * @brief Get field state
     * @return Non-const reference to state
     */
    FieldStateType &getState() { return *state_; }
    /**
     * @brief Get field state
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
    const bool is_subfield_; // Indicates whether this field is a sub-field of a
                             // rank > 0 tensor.
                             // Example: rank = 1 tensor
                             // tensor_component[0]: is_subfield_ = false
                             // tensor_component[1]: is_subfield_ = true
                             // tensor_component[n]: is_subfield_ = true
    FieldStateType *state_;  // Field state values/meta data

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

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr size_t Field<T, Entity, DIM, State, Alloc>::Rank;

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr size_t Field<T, Entity, DIM, State, Alloc>::NComponents;

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr Cubism::EntityType Field<T, Entity, DIM, State, Alloc>::EntityType;

template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr Cubism::FieldClass Field<T, Entity, DIM, State, Alloc>::Class;

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
    using MemoryOwner = typename BlockDataType::MemoryOwner;
    using ContainerType = std::vector<BaseType *>;

private:
    using FieldStateType = typename FieldType::FieldStateType;

public:
    /**
     * @brief Default constructor
     */
    FieldContainer() = default;

    /**
     * @brief Standard constructor for field container
     * @param n Number of components
     * @param r Index range spanned by field data
     * @param fs Field state initial value
     */
    FieldContainer(const size_t n,
                   const IndexRangeType &r,
                   const FieldStateType &fs = FieldStateType())
        : components_(n, nullptr)
    {
        for (size_t i = 0; i < n; ++i) {
            components_[i] = new FieldType(r, fs);
        }
    }

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
     * @brief Get container of raw data
     * @return STL vector with ``BaseType`` pointers
     */
    ContainerType &getContainer() { return components_; }

    /**
     * @brief Get container of raw data
     * @return ``const`` STL vector with ``BaseType`` pointers
     */
    const ContainerType &getContainer() const { return components_; }

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

/** @brief Generic tensor field
 * @tparam T Field data type
 * @tparam RANK Tensor rank
 * @tparam Entity Entity type
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator */
template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
class TensorField : public FieldContainer<Field<T, Entity, DIM, State, Alloc>>
{
public:
    using BaseType = FieldContainer<Field<T, Entity, DIM, State, Alloc>>;
    using typename BaseType::BlockDataType;
    using typename BaseType::DataType;
    using typename BaseType::FieldType; // scalar field sub-type
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using FieldStateType = typename FieldType::FieldStateType;
    using MemoryOwner = typename BlockDataType::MemoryOwner;

private:
    template <size_t B, size_t E>
    struct Power {
        static constexpr size_t value = B * Power<B, E - 1>::value;
    };
    template <size_t B>
    struct Power<B, 0> {
        static constexpr size_t value = 1;
    };
    template <typename U>
    using NestedVector = std::vector<std::vector<U>>; // simplifies construction
    using BaseType::components_;

public:
    static constexpr size_t Rank = RANK;
    static constexpr size_t NComponents =
        Power<IndexRangeType::Dim, RANK>::value;
    static constexpr Cubism::EntityType EntityType = BlockDataType::EntityType;
    static constexpr Cubism::FieldClass Class = Cubism::FieldClass::Tensor;
    static_assert(NComponents > 0, "Tensor has zero components");
    static_assert(IndexRangeType::Dim > 0, "DIM must be greater than zero");

    /** @brief Default constructor */
    TensorField() = delete;

    /**
     * @brief Main constructor to generate a tensor field
     * @param r Index range spanned by the field data
     * @param fs Field state initial value
     */
    explicit TensorField(const IndexRangeType &r,
                         const FieldStateType &fs = FieldStateType())
        : BaseType()
    {
        // first component owns the field state
        FieldType *first = new FieldType(r, fs);
        FieldStateType *pfs = &first->getState();
        components_.push_back(first);
        for (size_t i = 1; i < NComponents; ++i) {
            components_.push_back(new FieldType(r, pfs));
        }
    }

    /**
     * @brief Low-level field constructor (for higher rank tensors)
     * @param r Index range that spans the data
     * @param pfs Field state pointer (externally managed)
     */
    explicit TensorField(const IndexRangeType &r, FieldStateType *pfs)
        : BaseType()
    {
        for (size_t i = 0; i < NComponents; ++i) {
            components_.push_back(new FieldType(r, pfs));
        }
    }

    /**
     * @brief Low-level copy constructor for deep and shallow copies
     * @param tfc Tensor field container to be copied
     * @param o Memory ownership (``BlockDataType::MemoryOwner::Yes``: deep
     * copy)
     */
    TensorField(const TensorField &tfc, const MemoryOwner o) : BaseType()
    {
        // first component owns the field state
        FieldType *first = new FieldType(tfc[0], o);
        FieldStateType *pfs = &first->getState();
        components_.push_back(first);
        for (size_t i = 1; i < NComponents; ++i) {
            components_.push_back(new FieldType(tfc[i], o, pfs));
        }
    }

    /**
     * @brief Low-level copy constructor for deep and shallow copies (for higher
     * rank tensors)
     * @param tfc Tensor field container to be copied
     * @param o Memory ownership (``BlockDataType::MemoryOwner::Yes``: deep
     * copy)
     * @param pfs External field state
     */
    TensorField(const TensorField &tfc,
                const MemoryOwner o,
                const FieldStateType *pfs)
        : BaseType()
    {
        // first component owns the field state
        for (size_t i = 0; i < NComponents; ++i) {
            components_.push_back(new FieldType(tfc[i], o, pfs));
        }
    }

    /**
     * @brief Low-level constructor for external memory management
     * @param range_list Vector of index range
     * @param ptr_list Vector of block data pointer
     * @param bytes_list Number of bytes pointed to by pointer in ptr_list
     * @param state_list Vector of field state pointer
     * @param subfield Subfield component indicator
     *
     * The nested vector data structure is used to simplify the interface
     * between scalar fields, tensor fields and face field containers.  A
     * subfield component applies to higher rank tensors and face field
     * containers which share the field state among each others.
     */
    TensorField(const NestedVector<IndexRangeType> &range_list,
                const NestedVector<DataType *> &ptr_list,
                const NestedVector<size_t> &bytes_list,
                const NestedVector<FieldStateType *> &state_list,
                const bool subfield = false)
        : BaseType()
    {
        // outer
        assert(range_list.size() == 1);
        assert(ptr_list.size() == 1);
        assert(bytes_list.size() == 1);
        assert(state_list.size() == 1);
        // inner
        assert(range_list[0].size() == NComponents);
        assert(ptr_list[0].size() == NComponents);
        assert(bytes_list[0].size() == NComponents);
        assert(state_list[0].size() == NComponents);
        std::vector<bool> is_subfield(NComponents, true);
        if (!subfield) {
            is_subfield[0] = false;
        }
        for (size_t i = 0; i < NComponents; ++i) {
            const NestedVector<IndexRangeType> A(
                1,
                typename NestedVector<IndexRangeType>::value_type(
                    1, range_list[0][i]));
            const NestedVector<DataType *> B(
                1,
                typename NestedVector<DataType *>::value_type(1,
                                                              ptr_list[0][i]));
            const NestedVector<size_t> C(
                1,
                typename NestedVector<size_t>::value_type(1, bytes_list[0][i]));
            const NestedVector<FieldStateType *> D(
                1,
                typename NestedVector<FieldStateType *>::value_type(
                    1, state_list[0][i]));
            components_.push_back(new FieldType(A, B, C, D, is_subfield[i]));
        }
    }

    /**
     * @brief Copy constructor
     * @param c Tensor field to copy from
     */
    TensorField(const TensorField &c) : BaseType()
    {
        FieldType *first = new FieldType(c[0]);
        FieldStateType *pfs = &first->getState();
        components_.push_back(first);
        for (size_t i = 1; i < NComponents; ++i) {
            components_.push_back(
                new FieldType(c[i], first->getMemoryOwnership(), pfs));
        }
    }

    TensorField(TensorField &&c) = default;
    TensorField &operator=(const TensorField &rhs) = default;
    TensorField &operator=(TensorField &&c) = default;
    ~TensorField() = default;

    /**
     * @brief Get field state
     * @return Reference to state
     *
     * A tensor field shares one field state instance with all its components.
     */
    FieldStateType &getState() { return components_[0]->getState(); }

    /**
     * @brief Get field state
     * @return ``const`` reference to state
     *
     * A tensor field shares one field state instance with all its components.
     */
    const FieldStateType &getState() const
    {
        return components_[0]->getState();
    }

    /**
     * @brief Get memory ownership
     * @return Enumeration type describing the memory ownership
     */
    MemoryOwner getMemoryOwnership() const
    {
        return components_[0]->getMemoryOwnership();
    }

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

template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr size_t TensorField<T, RANK, Entity, DIM, State, Alloc>::Rank;

template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr size_t TensorField<T, RANK, Entity, DIM, State, Alloc>::NComponents;

template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr Cubism::EntityType
    TensorField<T, RANK, Entity, DIM, State, Alloc>::EntityType;

template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
constexpr Cubism::FieldClass
    TensorField<T, RANK, Entity, DIM, State, Alloc>::Class;

/** @brief Container class for all faces in a ``CUBISM_DIMENSION``-ional problem
 * @tparam TField Face field type (scalar or tensor)
 *
 * @rst
 * For ``CUBISM_DIMENSION`` in ``{1,2,3}``, the face field for faces with normal
 * in the ``X`` direction can be obtained with ``ff[0]`` or
 * ``ff[Cubism::Dir::X]`` for example, where ``ff`` is of type
 * ``FaceContainer``.
 * @endrst
 * */
template <typename TField>
class FaceContainer : public FieldContainer<TField>
{
public:
    using BaseType = FieldContainer<TField>;
    using FaceComponentType = TField; // scalar field or tensor field
    using typename BaseType::BlockDataType;
    using typename BaseType::DataType;
    using typename BaseType::FieldType; // scalar field sub-type
    using typename BaseType::IndexRangeType;
    using typename BaseType::MultiIndex;
    using FieldStateType = typename FieldType::FieldStateType;
    using MemoryOwner = typename BlockDataType::MemoryOwner;

private:
    template <typename U>
    using NestedVector = std::vector<std::vector<U>>; // simplifies construction
    using BaseType::components_;

public:
    static constexpr size_t Rank = FaceComponentType::Rank;
    static constexpr size_t NComponents = FaceComponentType::NComponents;
    static constexpr Cubism::EntityType EntityType = BlockDataType::EntityType;
    static constexpr Cubism::FieldClass Class =
        Cubism::FieldClass::FaceContainer;
    static_assert(
        BlockDataType::EntityType == Cubism::EntityType::Face,
        "FaceContainer: Entity type of field must be Cubism::EntityType::Face");
    static_assert(IndexRangeType::Dim > 0, "DIM must be greater than zero");

    /** @brief Default constructor */
    FaceContainer() = delete;

    /**
     * @brief Main constructor to generate a face container
     * @param cell_domain Index range spanned by the cell domain
     * @param fs Field state initial value
     */
    explicit FaceContainer(const IndexRangeType &cell_domain,
                           const FieldStateType &fs = FieldStateType())
        : BaseType()
    {
        const MultiIndex cells = cell_domain.getExtent(); // number of cells
        const IndexRangeType r0(cells + MultiIndex::getUnitVector(0));
        FaceComponentType *first = new FaceComponentType(r0, fs);
        FieldStateType *pfs = &first->getState();
        components_.push_back(first);
        for (size_t i = 1; i < IndexRangeType::Dim; ++i) {
            const IndexRangeType ri(cells + MultiIndex::getUnitVector(i));
            components_.push_back(new FaceComponentType(ri, pfs));
        }
    }

    /**
     * @brief Low-level copy constructor for deep and shallow copies
     * @param ffc Face field container to be copied
     * @param o Memory ownership (``BlockDataType::MemoryOwner::Yes``: deep
     * copy)
     */
    FaceContainer(const FaceContainer &ffc, const MemoryOwner o) : BaseType()
    {
        FaceComponentType *first = new FaceComponentType(ffc[0], o);
        FieldStateType *pfs = &first->getState();
        components_.push_back(first);
        for (size_t i = 1; i < IndexRangeType::Dim; ++i) {
            components_.push_back(new FaceComponentType(ffc[i], o, pfs));
        }
    }

    /**
     * @brief Low-level constructor for external memory management
     * @param range_list Vector of index range
     * @param ptr_list Vector of block data pointer
     * @param bytes_list Number of bytes pointed to by pointer in ptr_list
     * @param state_list Vector of field state pointer
     * @param subfield Subfield component indicator
     *
     * The nested vector data structure is used to simplify the interface
     * between scalar fields, tensor fields and face field containers.  A
     * subfield component applies to higher rank tensors and face field
     * containers which share the field state among each others.
     */
    FaceContainer(const NestedVector<IndexRangeType> &range_list,
                  const NestedVector<DataType *> &ptr_list,
                  const NestedVector<size_t> &bytes_list,
                  const NestedVector<FieldStateType *> &state_list,
                  const bool subfield = false)
        : BaseType()
    {
        // outer
        assert(range_list.size() == IndexRangeType::Dim);
        assert(ptr_list.size() == IndexRangeType::Dim);
        assert(bytes_list.size() == IndexRangeType::Dim);
        assert(state_list.size() == IndexRangeType::Dim);
        std::vector<bool> is_subfield(IndexRangeType::Dim, true);
        if (!subfield) {
            is_subfield[0] = false;
        }
        for (size_t i = 0; i < IndexRangeType::Dim; ++i) {
            // inner
            assert(range_list[i].size() == NComponents);
            assert(ptr_list[i].size() == NComponents);
            assert(bytes_list[i].size() == NComponents);
            assert(state_list[i].size() == NComponents);
            const NestedVector<IndexRangeType> A(
                1,
                typename NestedVector<IndexRangeType>::value_type(
                    range_list[i]));
            const NestedVector<DataType *> B(
                1, typename NestedVector<DataType *>::value_type(ptr_list[i]));
            const NestedVector<size_t> C(
                1, typename NestedVector<size_t>::value_type(bytes_list[i]));
            const NestedVector<FieldStateType *> D(
                1,
                typename NestedVector<FieldStateType *>::value_type(
                    state_list[i]));
            components_.push_back(
                new FaceComponentType(A, B, C, D, is_subfield[i]));
        }
    }

    /**
     * @brief Copy constructor
     * @param c Face container to copy from
     */
    FaceContainer(const FaceContainer &c) : BaseType()
    {
        FaceComponentType *first = new FaceComponentType(c[0]);
        FieldStateType *pfs = &first->getState();
        components_.push_back(first);
        for (size_t i = 1; i < IndexRangeType::Dim; ++i) {
            components_.push_back(
                new FaceComponentType(c[i], first->getMemoryOwnership(), pfs));
        }
    }

    FaceContainer(FaceContainer &&c) = default;
    FaceContainer &operator=(const FaceContainer &rhs) = default;
    FaceContainer &operator=(FaceContainer &&c) = default;
    ~FaceContainer() = default;

    /**
     * @brief Get field state
     * @return Reference to state
     *
     * A face field container shares one field state instance with all its
     * components.
     */
    FieldStateType &getState() { return components_[0]->getState(); }

    /**
     * @brief Get field state
     * @return ``const`` reference to state
     *
     * A face field container shares one field state instance with all its
     * components.
     */
    const FieldStateType &getState() const
    {
        return components_[0]->getState();
    }

    /**
     * @brief Get memory ownership
     * @return Enumeration type describing the memory ownership
     */
    MemoryOwner getMemoryOwnership() const
    {
        return components_[0]->getMemoryOwnership();
    }

    /**
     * @brief Face component access
     * @param i Direction index
     * @return Reference to face component
     */
    FaceComponentType &operator[](const size_t i)
    {
        assert(i < IndexRangeType::Dim);
        return *components_[i];
    }

    /**
     * @brief Face component access
     * @param i Direction index
     * @return ``const`` reference to face component
     */
    const FaceComponentType &operator[](const size_t i) const
    {
        assert(i < IndexRangeType::Dim);
        return *components_[i];
    }

    /**
     * @brief Face component access
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Direction of face
     * @return Reference to face component
     */
    template <typename Dir>
    FaceComponentType &operator[](const Dir &t)
    {
        return this->operator[](static_cast<size_t>(t));
    }

    /**
     * @brief Face component access
     * @tparam Dir Special type that defines a cast to ``size_t``
     * @param t Direction of face
     * @return ``const`` reference to face component
     */
    template <typename Dir>
    const FaceComponentType &operator[](const Dir &t) const
    {
        return this->operator[](static_cast<size_t>(t));
    }
};

template <typename TField>
constexpr size_t FaceContainer<TField>::Rank;

template <typename TField>
constexpr size_t FaceContainer<TField>::NComponents;

template <typename TField>
constexpr Cubism::EntityType FaceContainer<TField>::EntityType;

template <typename TField>
constexpr Cubism::FieldClass FaceContainer<TField>::Class;

/**
 * @brief Field view type
 * @tparam TField Field type
 *
 * @rst
 * Provides a view (shallow copy) for scalar fields, tensor fields or face field
 * containers.  The corresponding field interface is inherited.

 * .. note:: A view type never owns memory.

 * @endrst
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
    using MemoryOwner = typename BaseType::BlockDataType::MemoryOwner;

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

/**
 * @brief Basic cell-centered data field
 * @tparam T Data type (must be POD)
 * @tparam DIM Field dimension
 * @tparam State State type
 * @tparam Alloc Allocator type
 */
template <typename T,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using CellField = Field<T, EntityType::Cell, DIM, State, Alloc>;

/**
 * @brief Basic node-centered data field
 * @tparam T Data type (must be POD)
 * @tparam DIM Field dimension
 * @tparam State State type
 * @tparam Alloc Allocator type
 */
template <typename T,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using NodeField = Field<T, EntityType::Node, DIM, State, Alloc>;

/**
 * @brief Basic face-centered data field.
 * @tparam T Data type (must be POD)
 * @tparam DIM Field dimension
 * @tparam State State type
 * @tparam Alloc Allocator type
 *
 * @rst
 * Faces are stored individually for the dimensionality specified by
 * ``CUBISM_DIMENSION`` at compile time.  See the ``FaceContainer`` type for a
 * container of size ``CUBISM_DIMENSION``.
 * @endrst
 */
template <typename T,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using FaceField = Field<T, EntityType::Face, DIM, State, Alloc>;

/** @brief Convenience type for vector fields
 * @tparam T Field data type
 * @tparam Entity Entity type
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator */
template <typename T,
          Cubism::EntityType Entity,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
using VectorField = TensorField<T, 1, Entity, DIM, State, Alloc>;

/** @brief Field type factory for tensor fields
 * @tparam T Field data type
 * @tparam RANK Tensor rank
 * @tparam Entity Entity type
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator */
template <typename T,
          size_t RANK,
          Cubism::EntityType Entity,
          size_t DIM = CUBISM_DIMENSION,
          typename State = FieldState,
          template <typename> class Alloc = AlignedBlockAllocator>
struct FieldTypeFactory {
    using Type = TensorField<T, RANK, Entity, DIM, State, Alloc>;
};

/** @brief Field type factory for face tensor fields
 * @tparam T Field data type
 * @tparam RANK Tensor rank
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator */
template <typename T,
          size_t RANK,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
struct FieldTypeFactory<T, RANK, Cubism::EntityType::Face, DIM, State, Alloc> {
    using Type = FaceContainer<
        TensorField<T, RANK, Cubism::EntityType::Face, DIM, State, Alloc>>;
};

/** @brief Field type factory for scalar fields
 * @tparam T Field data type
 * @tparam Entity Entity type
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator */
template <typename T,
          Cubism::EntityType Entity,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
struct FieldTypeFactory<T, 0, Entity, DIM, State, Alloc> {
    using Type = Field<T, Entity, DIM, State, Alloc>;
};

/** @brief Field type factory for face scalar fields
 * @tparam T Field data type
 * @tparam DIM Field dimension
 * @tparam State Field state type
 * @tparam Alloc Memory allocator */
template <typename T,
          size_t DIM,
          typename State,
          template <typename>
          class Alloc>
struct FieldTypeFactory<T, 0, Cubism::EntityType::Face, DIM, State, Alloc> {
    using Type =
        FaceContainer<Field<T, Cubism::EntityType::Face, DIM, State, Alloc>>;
};

NAMESPACE_END(Block)
NAMESPACE_END(Cubism)

#endif /* FIELD_H_7NZ0QFMC */
