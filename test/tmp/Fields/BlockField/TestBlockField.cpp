// File       : TestBlockField.cpp
// Created    : Sat Apr 13 2019 10:45:37 PM (+0200)
// Author     : Fabian Wermelinger
// Description: BlockField test
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "TestBlockField.h"

using namespace std;
using namespace Cubism;

// Block size
#define BSX 32
#define BSY 64
#define BSZ 128

template <typename T>
using FieldCell = BlockField::FieldCell<T, BSX, BSY, BSZ>;

template <typename T>
using FieldNode = BlockField::FieldNode<T, BSX, BSY, BSZ>;

template <typename T>
using FieldFaceX = BlockField::FieldFaceX<T, BSX, BSY, BSZ>;

template <typename T>
using FieldFaceY = BlockField::FieldFaceY<T, BSX, BSY, BSZ>;

template <typename T>
using FieldFaceZ = BlockField::FieldFaceZ<T, BSX, BSY, BSZ>;

using FOwner = FieldCell<DataType>;
using FProxy = BlockField::FieldProxy<FOwner>;
using TReal = FProxy::DataType;

using FCell = FieldCell<DataType>;
using FNode = FieldNode<DataType>;

static void constructFields()
{
    tell("TEST: constructFields");
    Timer t0;
    ostringstream msg;

    t0.start();

    // owner
    FOwner fo;
    setFieldValue(fo, 1);
    assert(sumField(fo) == fo.getBlockSize());

    {
        // proxy
        // FProxy fp; // compile time error
        FProxy fp_0(fo.getBlockPtr()); // from block ptr
        assert(sumField(fp_0) == fp_0.getBlockSize());
    }
    assert(fo.getBlockPtr() != nullptr);

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}

static void copyConstructFields()
{
    tell("TEST: copyConstructFields");
    Timer t0;
    ostringstream msg;

    t0.start();

    // owner
    FOwner fo_0;
    setFieldValue(fo_0, 1);
    FOwner fo_1(fo_0);
    assert(fo_0.getBlockPtr() != fo_1.getBlockPtr());
    assert(sumField(fo_1) == sumField(fo_0));

    // proxy
    {
        FProxy fp_0(fo_0.getBlockPtr());                  // from block ptr
        FProxy fp_1(fo_0);                                // from field
        assert(fp_0.getBlockPtr() == fo_0.getBlockPtr());
        assert(fp_0.getBlockPtr() == fp_1.getBlockPtr());

        FProxy fp_2(fp_0);
        assert(fp_0.getBlockPtr() == fp_2.getBlockPtr());
    }
    assert(fo_0.getBlockPtr() != nullptr);

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}

static void moveConstructFields()
{
    tell("TEST: moveConstructFields");
    Timer t0;
    ostringstream msg;

    t0.start();

    // owner
    FOwner fo_0;
    setFieldValue(fo_0, 1);
    const void *pfo_0 = fo_0.getBlockPtr();

    FOwner fo_1(std::move(fo_0));
    assert(fo_1.getBlockPtr() == pfo_0);
    assert(sumField(fo_1) == fo_1.getBlockSize());

    {
        // proxy
        // FProxy fp_0(std::move(fo_1)); // compile time error
        // FProxy fp_1(std::move(fp_0)); // compile time error
    }

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}

static void assignFields()
{
    tell("TEST: assignFields");
    Timer t0;
    ostringstream msg;

    t0.start();

    // owner
    FOwner fo_0;
    setFieldValue(fo_0, 1);
    assert(sumField(fo_0) == fo_0.getBlockSize());
    FOwner fo_1;
    setFieldValue(fo_1, 2);
    assert(sumField(fo_1) == 2 * fo_1.getBlockSize());
    fo_1 = fo_0;
    assert(fo_0.getBlockPtr() != fo_1.getBlockPtr());
    assert(sumField(fo_1) == fo_1.getBlockSize());
    const double towner = t0.stop();

    t0.start();
    {
        // proxy
        FProxy fp_0(fo_0.getBlockPtr());
        FProxy fp_1(fo_1);
        setFieldValue(fp_1, 3);
        assert(sumField(fp_1) == 3 * fo_1.getBlockSize());
        fp_1 = fp_0; // proxies (operates on data, deep)
        assert(fp_0.getBlockPtr() != fp_1.getBlockPtr());
        assert(sumField(fp_1) == fp_1.getBlockSize());

        fp_0 = fo_1; // assign field to proxy (shallow)
        assert(fp_0.getBlockPtr() == fo_1.getBlockPtr());
        assert(sumField(fp_0) == fo_1.getBlockSize());
    }
    const double tproxy = t0.stop();
    assert(fo_0.getBlockPtr() != nullptr);
    assert(fo_1.getBlockPtr() != nullptr);

    msg << "Took: Owner fields = " << towner
        << " sec; Proxy fields = " << tproxy << " sec" << '\n';
    tell(msg.str());
}

static void moveAssignFields()
{
    tell("TEST: moveAssignFields");
    Timer t0;
    ostringstream msg;

    t0.start();

    // owner
    FOwner fo_0;
    setFieldValue(fo_0, 1);
    const void *pfo_0 = fo_0.getBlockPtr();

    FOwner fo_1;
    setFieldValue(fo_1, 2);
    fo_1 = std::move(fo_0);
    assert(fo_1.getBlockPtr() == pfo_0);
    assert(sumField(fo_1) == fo_1.getBlockSize());

    {
        // proxy
        // FProxy fp_0(
        //     fo_0.getBlockPtr()); // run time error (dereference a nullptr)
        FProxy fp_0(fo_1);
        FProxy fp_1(fo_1);
        // fp_1 = std::move(fp_0); // compile time error
    }
    assert(fo_0.getBlockPtr() == nullptr);
    assert(fo_1.getBlockPtr() != nullptr);

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}

static void basePointer()
{
    tell("TEST: basePointer");
    Timer t0;
    ostringstream msg;

    t0.start();

    FCell fc;
    FNode fn;
    setFieldValue(fc, 1);
    setFieldValue(fn, 1);
    assert(sumField(fc) == fc.getBlockSize());
    assert(sumField(fn) == fn.getBlockSize());

    // compile time error
    // typename FCell::BaseType *base_c = &fn;
    // typename FNode::BaseType *base_n = &fc;

    typename FCell::BaseType *base_c = &fc;
    typename FNode::BaseType *base_n = &fn;
    assert(base_c->getBlockSize() == fc.getBlockSize());
    assert(base_n->getBlockSize() == fn.getBlockSize());
    assert(base_c->getBlockPtr() == fc.getBlockPtr());
    assert(base_n->getBlockPtr() == fn.getBlockPtr());

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}

static void staticMember()
{
    tell("TEST: staticMember");
    Timer t0;
    ostringstream msg;

    t0.start();
    printStaticMember<FieldCell<DataType>>(msg);
    printStaticMember<FieldNode<DataType>>(msg);
    printStaticMember<FieldFaceX<DataType>>(msg);
    printStaticMember<FieldFaceY<DataType>>(msg);
    printStaticMember<FieldFaceZ<DataType>>(msg);
    printStaticMember<FProxy>(msg);
    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}

extern void blockIndexing();

int main(void)
{
    assert(sizeof(DataType) == sizeof(TReal));
    cout << "sizeof(TReal)    = " << sizeof(TReal) << '\n';
    cout << "sizeof(DataType) = " << sizeof(DataType) << '\n';
    cout << "sizeof(FOwner)   = " << sizeof(FOwner) << '\n';
    cout << "sizeof(FProxy)   = " << sizeof(FProxy) << '\n';
    cout << '\n';

    constructFields();
    copyConstructFields();
    moveConstructFields();
    assignFields();
    moveAssignFields();

    basePointer();

    staticMember();

    blockIndexing();

    return 0;
}
