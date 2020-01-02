// File       : Histogram.h
// Created    : Mon Dec 23 2019 10:57:51 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Histogram profiling on given MPI communicator
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef HISTOGRAM_H_6SPETKCI
#define HISTOGRAM_H_6SPETKCI

#include "Common.h"
#include "Util/Sampler.h"

#include <mpi.h>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

class Histogram : public Sampler
{
public:
    Histogram(const MPI_Comm comm,
              const std::string &name,
              const bool active = true)
        : Sampler(active), comm_(comm), name_(name)
    {
    }
    ~Histogram()
    {
        if (active_) {
            consolidate_();
        }
    }

    Histogram(const Histogram &c) = delete;
    Histogram &operator=(const Histogram &c) = delete;

private:
    const MPI_Comm comm_;
    const std::string name_;

    void consolidate_();
    void homogenizeCollection_();
};

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* HISTOGRAM_H_6SPETKCI */
