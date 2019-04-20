// File       : FieldIndexing.cpp
// Created    : Sun Apr 14 2019 03:26:12 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Indexing benchmark for field types
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "TestBlockField.h"

#include <functional>

using namespace std;
using namespace Cubism;

#define SAMPLES 10
#define FLOP 1

template <size_t BX, size_t BY, size_t BZ>
static void profileIndexing(ostream &s)
{
    using TestType = BlockField::FieldCell<DataType, BX, BY, BZ>;
    TestType block;
    setFieldValue(block, 1.0);

    auto linear = [&block](void) {
        DataType sum = 0.0;
        for (size_t i = 0; i < block.getBlockSize(); ++i) {
            sum += block[i]; // 1 Flop
        }
        return sum;
    };

    auto ijk = [&block](void) {
        DataType sum = 0.0;
        for (size_t iz = 0; iz < TestType::BlockDimZ; ++iz)
            for (size_t iy = 0; iy < TestType::BlockDimY; ++iy)
                for (size_t ix = 0; ix < TestType::BlockDimX; ++ix) {
                    sum += block(ix, iy, iz); // 1 Flop
                }
        return sum;
    };

    auto bench = [](std::function<DataType(void)> kern) {
        DataType res = 0.0;
        for (size_t i = 0; i < SAMPLES; ++i) {
            res += kern();
        }
        return res / SAMPLES;
    };

    Timer t0;

    s << "Size = " << BX << "x" << BY << "x" << BZ << " (" << BX * BY * BZ
      << ")" << '\n';

    { // linear indexing
        // warmup
        DataType res = linear();
        // profile
        t0.start();
        res += bench(linear);
        const double tlinear = t0.stop();
        s << "\tLinear indexing: result = " << static_cast<size_t>(res / 2)
          << "; took " << tlinear << " sec (" << SAMPLES << " samples)" << '\n';
        s << "\tAverage GFlop/s = "
          << (FLOP * SAMPLES * BX * BY * BZ) / tlinear * 1.0e-9 << '\n';
    }

    { // IJK indexing
        // warmup
        DataType res = ijk();
        // profile
        t0.start();
        res += bench(ijk);
        const double tijk = t0.stop();
        s << "\tIJK indexing: result = " << static_cast<size_t>(res / 2)
          << "; took " << tijk << " sec (" << SAMPLES << " samples)" << '\n';
        s << "\tAverage GFlop/s = "
          << (FLOP * SAMPLES * BX * BY * BZ) / tijk * 1.0e-9 << '\n';
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
    profileIndexing<512, 512, 512>(msg);

    msg << "Took: " << t0.stop() << " sec" << '\n';
    tell(msg.str());
}
