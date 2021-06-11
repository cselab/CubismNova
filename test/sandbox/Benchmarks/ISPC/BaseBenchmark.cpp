// File       : BaseBenchmark.cpp
// Created    : Wed Jun 09 2021 11:53:58 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Benchmark base class implementation
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "BaseBenchmark.h"

#include "Cubism/IO/Data.h"
#include "Cubism/Util/Timer.h"

#include <cmath>
#include <cstdio>

BaseBenchmark::BaseBenchmark(const int n_samples,
                             const int n_elements_per_dim,
                             const int stencil_start,
                             const int stencil_end,
                             const bool is_tensorial)
    : n_samples_(n_samples), field_(IndexRange(n_elements_per_dim))
{
    const Stencil stencil(stencil_start, stencil_end, is_tensorial);
    lab_.allocate(stencil, field_.getIndexRange());
}

void BaseBenchmark::writeTestData(const std::string label)
{
    Cubism::IO::DataWriteUniformHDF<Real>(
        label + "_field",
        "data",
        static_cast<typename Field::BlockDataType>(field_));
    Cubism::IO::DataWriteUniformHDF<Real>(
        label + "_lab",
        "data",
        static_cast<typename FieldLab::BlockDataType>(lab_));
}

void BaseBenchmark::writeInitTestData()
{
    init_();
    this->writeTestData("initial");
}

typename BaseBenchmark::Result
BaseBenchmark::benchmark_(const std::string &tag,
                          BaseBenchmark::Kernel kernel,
                          const int loop_flop,
                          const typename BaseBenchmark::Result *gold)
{
    // prepare the test data
    init_();

    // extract kernel arguments
    const auto extent = field_.getIndexRange().getExtent();
    Real *dst = field_.getData();
    const int x_pitch_dst = extent[0];
    const int xy_pitch_dst = extent[0] * extent[1];

    const auto lab_extent = lab_.getIndexRange().getExtent();
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
           xy_pitch_dst,
           char_spacing_);

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
               xy_pitch_dst,
               char_spacing_);
        samples[i] = t.stop();
    }

    const int total_elements = field_.size();
    return report_(tag, total_elements * loop_flop, samples, gold);
}

typename BaseBenchmark::Result
BaseBenchmark::report_(const std::string &tag,
                       const int flop,
                       const std::vector<double> &samples,
                       const typename BaseBenchmark::Result *gold)
{
    Result res;
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

    std::printf("\033[1;35m%-16s\033[0m avg:%.4ems std:%.4ems min:%.4ems "
                "max:%.4ems samples:%d\n",
                tag.c_str(),
                1.0e3 * res.mean,
                1.0e3 * res.sdev,
                1.0e3 * res.min,
                1.0e3 * res.max,
                n_samples);
    std::printf(
        "%-16s \033[1;34mGflop/s:%.4e\033[0m", "", flop / res.mean * 1.0e-9);
    if (gold != nullptr) {
        const double diff = res.mean - gold->mean;
        const double speedup = gold->mean / res.mean;
        std::printf(" [ref:%s diff:%.4ems \033[1;33mspeedup:%.3f\033[0m]\n",
                    (gold->tag).c_str(),
                    1.0e3 * diff,
                    speedup);
    } else {
        std::printf("\n");
    }

    return res;
}

void BaseBenchmark::init_()
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
    hinv_ = 1.0 / static_cast<Real>(extent[0]);

    // load the lab with periodic boundary
    auto periodic = [&](const MIndex &) -> Field & { return field_; };
    lab_.loadData(MIndex(0), periodic);
}
