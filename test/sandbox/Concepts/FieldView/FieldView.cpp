// File       : field_view.cpp
// Created    : Thu Apr 11 2019 03:47:19 PM (+0200)
// Author     : Fabian Wermelinger
// Description: FieldView: Type-safe collection of references to objects
//              (fields) with possibly different data type.  Generic approach is
//              preferred as many different types of fields may exist and
//              constraining them to use the same base class interface is too
//              rigid.
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include <cassert>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

#define SIZE 100

/// @brief Minimalistic representation of a field with specific type
///
/// @tparam T
template <typename T>
class Field
{
public:
    Field(std::string name) : v_(SIZE), name_(std::move(name))
    {
        std::cout << "Constructor Field<" << typeid(T).name() << ">(" << name_
                  << "): callid: " << ++ctor_count << '\n';
    }

    // disable copying
    Field(Field &&c) = delete;
    Field(const Field &c) = delete;
    Field &operator=(const Field &c) = delete;

    using TData = T;

    std::string name() const { return name_; }

    const T &operator[](const size_t i) const
    {
        assert(i < SIZE);
        return v_[i];
    }

    T &operator[](const size_t i)
    {
        assert(i < SIZE);
        return v_[i];
    }

private:
    static size_t ctor_count;
    std::vector<T> v_;
    std::string name_;
};

template <typename T>
size_t Field<T>::ctor_count = 0;

/// @brief FieldView prototype
///
/// @tparam TFields
template <typename... TFields>
class FieldView
{
public:
    FieldView(TFields &... fields) // lvalues only
        : fields_(std::tuple<TFields &&...>(std::forward<TFields>(fields)...))
    {
    }

    // disable copying
    FieldView(FieldView &&c) = delete;
    FieldView(const FieldView &c) = delete;
    FieldView &operator=(const FieldView &c) = delete;

    // XXX: [fabianw@mavt.ethz.ch; 2019-04-11] This iterator applies to fields
    // which hold a collection of 'blocks'.  In this case the iterator shall
    // return a reference to the block pointed to by the individual field
    // iterators.  Alternatively, a utility function outside FieldView might be
    // a better solution: FieldUtils::zip(FieldA fa, FieldB fb, ...)
    // class iterator
    // {
    // public:
    //     iterator() {}
    //     ~iterator() {}

    // private:
    // };

    template <size_t id, typename FieldType>
    FieldType &getField()
    {
        static_assert(id < sizeof...(TFields), "Field id out of bounds");
        return std::get<id>(fields_);
    }

    template <size_t id>
    std::string getFieldName()
    {
        static_assert(id < sizeof...(TFields), "Field id out of bounds");
        return std::get<id>(fields_).name();
    }

private:
    std::tuple<TFields &&...> fields_;
};

int main()
{
    struct AoS {
        double d0;
        int i0;
    };
    using FD = Field<double>;
    using FI = Field<int>;
    using FC = Field<char>;
    using FAoS = Field<AoS>;         // Array of structures layout
    using FVec2 = FieldView<FD, FD>; // combine 2 scalar fields

    // allocate some fields
    FD fd0("d0"), fd1("d1");
    FI fi0("i0");
    FC fc0("c0");
    FAoS faos("AoS");
    FVec2 fvd(fd0, fd1); // a reference to fd0 and fd1

    // initialize data
    for (size_t i = 0; i < SIZE; ++i) {
        fd0[i] = fd1[i] = fi0[i] = fc0[i] = i;
        faos[i].d0 = i;
        faos[i].i0 = -static_cast<int>(i);
    }

    // container of various types (just lvalues here)
    FieldView<FD, FD, FI, FC, FVec2, FAoS> fview(fd0, fd1, fi0, fc0, fvd, faos);

    std::cout << "Field name: fview<0> = " << fview.template getFieldName<0>()
              << '\n';
    std::cout << "Field name: fview<1> = " << fview.template getFieldName<1>()
              << '\n';
    std::cout << "Field name: fview<2> = " << fview.template getFieldName<2>()
              << '\n';
    std::cout << "Field name: fview<3> = " << fview.template getFieldName<3>()
              << '\n';
    std::cout << "Field name: fview<4>.comp<0> = "
              << fview.template getField<4, FVec2>().template getFieldName<0>()
              << '\n';
    std::cout << "Field name: fview<4>.comp<1> = "
              << fview.template getField<4, FVec2>().template getFieldName<1>()
              << '\n';
    std::cout << "Field name: fview<5> = " << fview.template getFieldName<5>()
              << '\n';

    // FD &d0 = fview.template getField<0, FD>();
    // FD &d1 = fview.template getField<1, FD>();
    // FI &i0 = fview.template getField<2, FI>();
    // FC &c0 = fview.template getField<3, FC>();
    // FVec2 &fv = fview.template getField<4, FVec2>();
    // FAoS &fa = fview.template getField<5, FAoS>();
    auto &d0 = fview.template getField<0, FD>();
    auto &d1 = fview.template getField<1, FD>();
    auto &i0 = fview.template getField<2, FI>();
    auto &c0 = fview.template getField<3, FC>();
    auto &fv = fview.template getField<4, FVec2>();
    auto &fa = fview.template getField<5, FAoS>();

    // compile time error
    // auto &dummy = fview.template getField<100, FVec2>();

    std::cout << "lvalue (const): d0[0]  = " << d0[0] << '\n';
    std::cout << "lvalue (const): d1[1]  = " << d1[1] << '\n';
    std::cout << "lvalue (const): i0[2]  = " << i0[2] << '\n';
    std::cout << "lvalue (const): c0[97] = " << c0[97] << '\n';
    std::cout << "fv.comp0[0] = " << (fv.template getField<0, FD>())[0] << '\n';
    std::cout << "fv.comp1[1] = " << (fv.template getField<1, FD>())[1] << '\n';

    std::cout << "lvalue (assign 2):   d0[0]  = " << (d0[0] = 2) << '\n';
    std::cout << "lvalue (assign 3):   d1[1]  = " << (d1[1] = 3) << '\n';
    std::cout << "lvalue (assign 4):   i0[2]  = " << (i0[2] = 4) << '\n';
    std::cout << "lvalue (assign 'b'): c0[97] = " << (c0[97] = 'b') << '\n';

    std::cout << "fv.comp0[0] = " << (fv.template getField<0, FD>())[0] << '\n';
    std::cout << "fv.comp1[1] = " << (fv.template getField<1, FD>())[1] << '\n';

    std::cout << "fa[1].d0 = " << fa[1].d0 << '\n';
    std::cout << "fa[1].i0 = " << fa[1].i0 << '\n';

    return 0;
}
