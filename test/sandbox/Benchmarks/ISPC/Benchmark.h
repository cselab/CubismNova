// File       : Benchmark.h
// Created    : Wed Jun 09 2021 11:53:58 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Benchmark base class
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef BENCHMARK_H_LYAQO1H7
#define BENCHMARK_H_LYAQO1H7

#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldLab.h"
#include "Cubism/IO/Data.h"
#include "Cubism/Util/Timer.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

template <typename Real>
class Benchmark
{
public:

    using Field = Cubism::Block::CellField<Real>;
    using FieldLab = Cubism::Block::FieldLab<Field>;
    using IndexRange = typename Field::IndexRangeType;
    using Stencil = typename FieldLab::StencilType;
    // low-level signature of compute kernel
    using Kernel = void (*)(
        const int Nx,            // number of x-elements
        const int Ny,            // number of y-elements
        const int Nz,            // number of z-elements
        const Real *src,         // pointer to source data
        const int x_pitch_src,   // x-pitch in src memory layout
        const int xy_pitch_src,  // xy-pitch (slice) in src memory layout
        Real *dst,               // pointer to destination data
        const int x_pitch_dst,   // x-pitch in dst memory layout
        const int xy_pitch_dst); // xy-pitch (slice) in dst memory layout

    Benchmark() = delete;
    Benchmark(const int n_samples,
              const int n_elements_per_dim,
              const int stencil_start,
              const int stencil_end,
              const bool is_tensorial = false)
        : n_samples_(n_samples), field_(IndexRange(n_elements_per_dim))
    {
        const Stencil stencil(stencil_start, stencil_end, is_tensorial);
        lab_.allocate(stencil, field_.getIndexRange());
    }
    virtual ~Benchmark() {}

    // main interface
    virtual void run() {}

    void writeTestData()
    {
        init_();
        Cubism::IO::DataWriteUniformHDF<Real>(
            "field",
            "data",
            static_cast<typename Field::BlockDataType>(field_));
        Cubism::IO::DataWriteUniformHDF<Real>(
            "lab", "data", static_cast<typename FieldLab::BlockDataType>(lab_));
    }

protected:
    struct Result {
        double mean, sdev, min, max;
        std::string tag;
    };

    Result benchmark_(const std::string &tag,
                      Kernel kernel,
                      const double flop,
                      const Result *gold = nullptr)
    {
        // prepare the test data
        init_();

        // extract kernel arguments
        const auto extent = field_.getIndexRange().getExtent();
        Real *dst = field_.getData();
        const int x_pitch_dst = extent[0];
        const int xy_pitch_dst = extent[0] * extent[1];

        const auto lab_extent = lab_.getMaximumRange().getExtent();
        const Real *src = lab_.getInnerData();
        const int x_pitch_src = lab_extent[0];
        const int xy_pitch_src = lab_extent[0] * lab_extent[1];

        // warm-up
        kernel(extent[0],
               extent[1],
               extent[2],
               src,
               x_pitch_src,
               xy_pitch_src,
               dst,
               x_pitch_dst,
               xy_pitch_dst);

        // collect samples
        Cubism::Util::Timer t;
        std::vector<double> samples(n_samples_, 0.0);
        for (size_t i = 0; i < samples.size(); ++i) {
            t.start();
            kernel(extent[0],
                   extent[1],
                   extent[2],
                   src,
                   x_pitch_src,
                   xy_pitch_src,
                   dst,
                   x_pitch_dst,
                   xy_pitch_dst);
            samples[i] = t.stop();
        }

        return report_(tag, flop, samples, gold);
    }

private:
    const int n_samples_;
    Field field_;
    FieldLab lab_;

    Result report_(const std::string &tag,
                   const double flop,
                   const std::vector<double> &samples,
                   const Result *gold = nullptr)
    {
        Result res{0};
        res.tag = tag;
        const int n_samples = samples.size();

        res.mean = 0.0;
        res.min = HUGE_VAL;
        res.max = 0.0;
        for (double s : samples) {
            res.mean += s;
            res.min = (s < res.min) ? s : res.min;
            res.max = (s > res.max) ? s : res.max;
        }
        res.mean /= n_samples;

        res.sdev = 0.0;
        for (double s : samples) {
            res.sdev += (s - res.mean) * (s - res.mean);
        }
        res.sdev = std::sqrt(res.sdev / (n_samples - 1.0));

        std::printf("%-10s: avg:%.4e std:%.4e min:%.4e max:%.4e samples:%d\n",
                    tag.c_str(),
                    res.mean,
                    res.sdev,
                    res.min,
                    res.max,
                    n_samples);
        std::printf("%-10s: Gflop/s:%.4e", "", flop / res.mean * 1.0e-9);
        if (gold != nullptr) {
            const double rel = (res.mean - gold->mean) / gold->mean;
            const double sup = gold->mean / res.mean;
            std::printf("[%-10s: rel:%.3f%% speedup:%.3f]\n",
                        (gold->tag).c_str(),
                        rel,
                        sup);
        } else {
            std::printf("\n");
        }

        return res;
    }

    void init_()
    {
        using MIndex = typename Field::MultiIndex;

        // initialize field data
        const auto extent = field_.getIndexRange().getExtent();
        for (auto i : field_.getIndexRange()) {
            const Real x = (static_cast<Real>(i[0]) + 0.5) / extent[0];
            const Real y = (static_cast<Real>(i[1]) + 0.5) / extent[1];
            const Real z = (static_cast<Real>(i[2]) + 0.5) / extent[2];
            field_[i] = std::sin(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y) *
                        std::sin(2.0 * M_PI * z);
        }

        // load the lab with periodic boundary
        auto periodic = [&](const MIndex &) -> Field & { return field_; };
        lab_.loadData(MIndex(0), periodic);
    }
};

#endif /* BENCHMARK_H_LYAQO1H7 */
