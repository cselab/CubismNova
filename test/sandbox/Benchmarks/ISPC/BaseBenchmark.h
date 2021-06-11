// File       : BaseBenchmark.h
// Created    : Wed Jun 09 2021 11:53:58 AM (+0200)
// Author     : Fabian Wermelinger
// Description: Benchmark base class
// Copyright 2021 ETH Zurich. All Rights Reserved.
#ifndef BASEBENCHMARK_H_AGNZJOFH
#define BASEBENCHMARK_H_AGNZJOFH

#include "Common.h"

#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldLab.h"
#include "Cubism/Core/Vector.h"

#include <string>
#include <vector>

class BaseBenchmark
{
public:
    using RealVec = Cubism::Core::Vector<Real>;
    using Field = Cubism::Block::CellField<Real>;
    using FieldLab = Cubism::Block::FieldLab<Field>;
    using IndexRange = typename Field::IndexRangeType;
    using Stencil = typename FieldLab::StencilType;
    // low-level signature of compute kernel
    using Kernel = void (*)(
        const int Nx,           // number of x-elements
        const int Ny,           // number of y-elements
        const int Nz,           // number of z-elements
        const Real *src,        // pointer to source data
        const int x_pitch_src,  // x-pitch in src memory layout
        const int xy_pitch_src, // xy-pitch (slice) in src memory layout
        Real *dst,              // pointer to destination data
        const int x_pitch_dst,  // x-pitch in dst memory layout
        const int xy_pitch_dst, // xy-pitch (slice) in dst memory layout
        const Real factor);     // multiplication factor (optional)

    BaseBenchmark() = delete;
    BaseBenchmark(const int n_samples,
                  const int n_elements_per_dim,
                  const int stencil_start,
                  const int stencil_end,
                  const bool is_tensorial = false);
    virtual ~BaseBenchmark() {}

    // main interface
    virtual void run() {}
    void writeInitTestData();
    void writeTestData(const std::string label = "");

protected:
    const int n_samples_;
    Field field_;
    Field gold_values_;
    FieldLab lab_;
    Real char_spacing_; // characteristic grid spacing

    struct Result {
        double mean, sdev, min, max, Gflop_per_second;
        std::string tag;
    };

    struct Error {
        double L1, L2, Linf;
    };

    virtual Error applicationSpecificError_(const Result *)
    {
        Error err;
        err.L1 = 0.0;
        err.L2 = 0.0;
        err.Linf = 0.0;
        return err;
    }
    virtual void benchmarkCustom_(const Result * = nullptr) {}
    Result benchmark_(const std::string &tag,
                      Kernel kernel,
                      const int loop_flop,
                      const Result *gold = nullptr);
    Error error_(const Result *gold)
    {
        if (gold == nullptr) {
            // This call initiates the reference.  Save the computed values for
            // relative error comparison due to different instructions and
            // instruction order of the benchmark kernels
            gold_values_ = field_;
        }
        return applicationSpecificError_(gold);
    }
    Result report_(const std::string &tag,
                   const int flop,
                   const std::vector<double> &samples,
                   const Error error,
                   const Result *gold = nullptr);
    void init_();
    Real f_(const Real x, const Real y, const Real z);
    RealVec gradf_(const Real x, const Real y, const Real z);
    Real divgradf_(const Real x, const Real y, const Real z);
};

#endif /* BASEBENCHMARK_H_AGNZJOFH */
