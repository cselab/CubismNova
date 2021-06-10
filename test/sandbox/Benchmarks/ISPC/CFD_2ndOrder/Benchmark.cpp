// File       : Benchmark.h
// Created    : Thu Jun 10 2021 09:41:23 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Centered finite differences 2nd order
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "CFD_2ndOrder/Benchmark.h"
#include "CFD_2ndOrder/Kernels.h"

#include "Cubism/Util/Timer.h"
#include <cstdio>

namespace CFD
{
namespace Order2
{
void Benchmark::run()
{
    printf("RUNNING: CFD_2ndOrder (2nd derivative)\n");

    Result gold = this->benchmark_("GOLD", Gold::ddx, Gold::loop_flop_ddx);

    this->benchmark_(
        "AUTOVEC_GOLD", Gold::ddxTreeVec, Gold::loop_flop_ddx_tree_vec, &gold);

    this->benchmarkCustom_(&gold);

    this->benchmark_(
        "ISPC_SSE2", ispc::ddxISPCsse2, ispc::loop_flop_ddx, &gold);

    this->benchmark_(
        "ISPC_SSE4", ispc::ddxISPCsse4, ispc::loop_flop_ddx, &gold);

    this->benchmark_("ISPC_AVX", ispc::ddxISPCavx, ispc::loop_flop_ddx, &gold);

    this->benchmark_(
        "ISPC_AVX2", ispc::ddxISPCavx2, ispc::loop_flop_ddx, &gold);
}

void Benchmark::benchmarkCustom_(const Result *gold)
{
    // prepare the test data
    init_();
    const int total_elements = field_.size();

    { // naive multi-index implementation
        // warm-up
        naiveMultiIndex_();

        // collect samples
        Cubism::Util::Timer t;
        std::vector<double> samples(n_samples_, 0.0);
        for (size_t i = 0; i < samples.size(); ++i) {
            t.start();
            naiveMultiIndex_();
            samples[i] = t.stop();
        }
        report_("NAIVE_MIDX", total_elements * naive_flop_, samples, gold);
    }

    { // compiler auto-vectorized multi-index implementation
        // warm-up
        naiveMultiIndexTreeVec_();

        // collect samples
        Cubism::Util::Timer t;
        std::vector<double> samples(n_samples_, 0.0);
        for (size_t i = 0; i < samples.size(); ++i) {
            t.start();
            naiveMultiIndexTreeVec_();
            samples[i] = t.stop();
        }
        report_("AUTOVEC_MIDX", total_elements * naive_flop_, samples, gold);
    }

    { // naive explicit 3-loop implementation
        // warm-up
        naive3DIndex_();

        // collect samples
        Cubism::Util::Timer t;
        std::vector<double> samples(n_samples_, 0.0);
        for (size_t i = 0; i < samples.size(); ++i) {
            t.start();
            naive3DIndex_();
            samples[i] = t.stop();
        }
        report_("NAIVE_3DIDX", total_elements * naive_flop_, samples, gold);
    }

    { // compiler auto-vectorized explicit 3-loop implementation
        // warm-up
        naive3DIndexTreeVec_();

        // collect samples
        Cubism::Util::Timer t;
        std::vector<double> samples(n_samples_, 0.0);
        for (size_t i = 0; i < samples.size(); ++i) {
            t.start();
            naive3DIndexTreeVec_();
            samples[i] = t.stop();
        }
        report_("AUTOVEC_3DIDX", total_elements * naive_flop_, samples, gold);
    }
}

#define NAIVE_MULTI_INDEX()                                                    \
    do {                                                                       \
        const MIndex ix{1, 0, 0};                                              \
        const MIndex iy{0, 1, 0};                                              \
        const MIndex iz{0, 0, 1};                                              \
        const Real fac = hinv_ * hinv_;                                        \
        for (const auto &i : field_.getIndexRange()) {                         \
            field_[i] = fac * (lab_[i - ix] + lab_[i + ix] + lab_[i - iy] +    \
                               lab_[i + iy] + lab_[i - iz] + lab_[i + iz] -    \
                               static_cast<Real>(6) * lab_[i]);                \
        }                                                                      \
    } while (0)

#if defined(__GNUC__) || defined(__clang__)
// no auto vectorization
__attribute__((optimize("no-tree-vectorize")))
#endif
void Benchmark::naiveMultiIndex_()
{
    using MIndex = typename IndexRange::MultiIndex;
    NAIVE_MULTI_INDEX();
}

void Benchmark::naiveMultiIndexTreeVec_()
{
    using MIndex = typename IndexRange::MultiIndex;
    NAIVE_MULTI_INDEX();
}

#undef NAIVE_MULTI_INDEX

#define NAIVE_3D_INDEX()                                                       \
    do {                                                                       \
        const Real fac = hinv_ * hinv_;                                        \
        const MIndex extent = field_.getIndexRange().getExtent();              \
        for (Index iz = 0; iz < extent[2]; ++iz) {                             \
            for (Index iy = 0; iy < extent[1]; ++iy) {                         \
                for (Index ix = 0; ix < extent[0]; ++ix) {                     \
                    field_(ix, iy, iz) =                                       \
                        fac * (lab_(ix - 1, iy, iz) + lab_(ix + 1, iy, iz) +   \
                               lab_(ix, iy - 1, iz) + lab_(ix, iy + 1, iz) +   \
                               lab_(ix, iy, iz - 1) + lab_(ix, iy, iz + 1) -   \
                               static_cast<Real>(6) * lab_(ix, iy, iz));       \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

#if defined(__GNUC__) || defined(__clang__)
// no auto vectorization
__attribute__((optimize("no-tree-vectorize")))
#endif
void Benchmark::naive3DIndex_()
{
    using MIndex = typename IndexRange::MultiIndex;
    using Index = typename MIndex::DataType;
    NAIVE_3D_INDEX();
}

void Benchmark::naive3DIndexTreeVec_()
{
    using MIndex = typename IndexRange::MultiIndex;
    using Index = typename MIndex::DataType;
    NAIVE_3D_INDEX();
}

#undef NAIVE_3D_INDEX

} // namespace Order2
} // namespace CFD
