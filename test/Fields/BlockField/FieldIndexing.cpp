// File       : FieldIndexing.cpp
// Created    : Sun Apr 14 2019 03:26:12 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Indexing benchmark for field types
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "TestBlockField.h"

#include <functional>

using namespace std;
using namespace Cubism;

template <size_t BX, size_t BY, size_t BZ>
static void profileIndexing(ostream &s)
{
    using TestType = BlockField::FieldCell<DataType, BX, BY, BZ>;
    TestType block;
    setFieldValue(block, 1.0);

    auto linear = [&block](void) {
        double sum = 0.0;
        for (size_t i = 0; i < block.getBlockSize(); ++i) {
            sum += block[i];
        }
        return sum;
    };

    auto ijk = [&block](void) {
        double sum = 0.0;
        for (size_t iz = 0; iz < TestType::BlockDimZ; ++iz)
            for (size_t iy = 0; iy < TestType::BlockDimY; ++iy)
                for (size_t ix = 0; ix < TestType::BlockDimX; ++ix) {
                    sum += block(ix, iy, iz);
                }
        return sum;
    };

    auto bench = [](std::function<double(void)> kern) {
        double res = 0.0;
        for (size_t i = 0; i < 10; ++i) {
            res += kern();
        }
        return res / 10;
    };

    Timer t0;

    s << "Size = " << BX << "x" << BY << "x" << BZ << " (" << BX * BY * BZ
      << ")" << '\n';

    {
        // warmup
        double res = linear();
        // profile
        t0.start();
        res += bench(linear);
        const double tlinear = t0.stop();
        s << "\tLinear indexing: result = " << static_cast<size_t>(res / 2)
          << "; took " << tlinear << " sec (10 samples)" << '\n';
    }

    {
        // warmup
        double res = ijk();
        // profile
        t0.start();
        res += bench(ijk);
        const double tijk = t0.stop();
        s << "\tIJK indexing: result = " << static_cast<size_t>(res / 2)
          << "; took " << tijk << " sec (10 samples)" << '\n';
    }
}

void blockIndexing()
{
    tell("TEST: blockIndexing");
    Timer t0;
    ostringstream msg;

    t0.start();
    profileIndexing<4, 4, 4>(msg);
    profileIndexing<16, 16, 16>(msg);
    profileIndexing<32, 32, 32>(msg);
    profileIndexing<64, 128, 256>(msg);
    profileIndexing<1024, 1024, 1024>(msg);

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}
