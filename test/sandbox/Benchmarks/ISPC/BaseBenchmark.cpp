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
    : n_samples_(n_samples), field_(IndexRange(n_elements_per_dim)),
      gold_values_(IndexRange(n_elements_per_dim))
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

BaseBenchmark::Result
BaseBenchmark::benchmark_(const std::string &tag,
                          BaseBenchmark::Kernel kernel,
                          const int loop_flop,
                          const BaseBenchmark::Result *gold)
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

    // compute error
    const Error error = error_(gold);

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
    return report_(tag, total_elements * loop_flop, samples, error, gold);
}

BaseBenchmark::Result BaseBenchmark::report_(const std::string &tag,
                                             const int flop,
                                             const std::vector<double> &samples,
                                             const BaseBenchmark::Error error,
                                             const BaseBenchmark::Result *gold)
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
    std::string error_string;
    if (gold == nullptr) {
        error_string = "Error (exact solution)";
    } else {
        error_string = "Error (relative to \033[1;35m" + gold->tag + "\033[0m)";
    }
    // std::printf("%-16s \033[1;31mError:\033[0m L1=%.6e L2=%.6e Linf=%.6e\n",
    std::printf("%-16s %s: L1=%.6e L2=%.6e Linf=%.6e\n",
                "",
                error_string.c_str(),
                error.L1,
                error.L2,
                error.Linf);
    res.Gflop_per_second = flop / res.mean * 1.0e-9;
    std::printf(
        "%-16s \033[1;34mGflop/s:%.4e\033[0m", "", res.Gflop_per_second);
    if (gold != nullptr) {
        const double diff = res.mean - gold->mean;
        const double speedup = res.Gflop_per_second / gold->Gflop_per_second;
        std::printf(" [ref:\033[1;35m%s\033[0m diff:%.4ems "
                    "\033[1;33mspeedup:%.3f\033[0m]\n",
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
    const MIndex extent = field_.getIndexRange().getExtent();
    char_spacing_ = 1.0 / static_cast<Real>(extent[0]);
    for (auto i : field_.getIndexRange()) {
        const Real x = char_spacing_ * (static_cast<Real>(i[0]) + 0.5);
        const Real y = char_spacing_ * (static_cast<Real>(i[1]) + 0.5);
        const Real z = char_spacing_ * (static_cast<Real>(i[2]) + 0.5);
        field_[i] = f_(x, y, z);
    }

    // load the lab with periodic boundary
    auto periodic = [&](const MIndex &) -> Field & { return field_; };
    lab_.loadData(MIndex(0), periodic);
}

Real BaseBenchmark::f_(const Real x, const Real y, const Real z)
{
    const Real fac = 2.0 * M_PI;
    return std::sin(fac * x) * std::cos(fac * y) * std::sin(fac * z);
}

BaseBenchmark::RealVec
BaseBenchmark::gradf_(const Real x, const Real y, const Real z)
{
    const Real fac = 2.0 * M_PI;
    return RealVec{
        fac * std::cos(fac * x) * std::cos(fac * y) * std::sin(fac * z),
        -fac * std::sin(fac * x) * std::sin(fac * y) * std::sin(fac * z),
        fac * std::sin(fac * x) * std::cos(fac * y) * std::cos(fac * z)};
}

Real BaseBenchmark::divgradf_(const Real x, const Real y, const Real z)
{
    const Real fac = -12.0 * M_PI * M_PI;
    return fac * f_(x, y, z);
}
