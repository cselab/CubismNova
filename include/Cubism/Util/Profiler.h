// File       : Profiler.h
// Created    : Mon Dec 23 2019 12:27:37 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Profiling agent
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef PROFILER_H_LQ95Z7KM
#define PROFILER_H_LQ95Z7KM

#include "Cubism/Common.h"
#include "Cubism/Util/Sampler.h"
#include <cassert>
#include <map>
#include <mpi.h>
#include <stack>
#include <string>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

/**
 * @ingroup Util MPI
 * @brief Runtime profiler
 *
 * @rst
 * Used to collect runtime samples for a code section that is enclosed
 * by the ``push()`` and ``pop()`` methods.  Used for profiling.
 * @endrst
 * */
class Profiler : public Sampler
{
    struct Accumulator {
        size_t total_samples;
        size_t batch_samples;
        double total_time_mean;
        double total_time_accu;
    };

public:
    /**
     * @brief Main constructor
     * @param name Name of the profiler
     * @param comm MPI communicator associated to the profiler
     */
    Profiler(const std::string &name = "default",
             const MPI_Comm comm = MPI_COMM_WORLD)
        : Sampler(), comm_(comm), name_(name), batch_count_(0)
    {
    }

    Profiler(const Profiler &c) = delete;
    Profiler &operator=(const Profiler &c) = delete;

    /**
     * @brief Push a new profiling agent on the stack
     * @param name Name of the agent
     *
     * Starts a timer for a new profiling agent.
     */
    void push(const std::string &name)
    {
        agents_.push(name);
        const auto it = agents_all_.find(name);
        if (it == agents_all_.end()) {
            agents_all_[name] = {};
        }
        this->seedSample();
    }

    /** @brief Pop the top profiling agent
     *
     * Pops the top (most recent) profiling agent from the stack and collects
     * the time measurement.
     * */
    void pop()
    {
        assert(!agents_.empty());
        this->collectSample(agents_.top());
        agents_.pop();
    }

    /** @brief Print profiling report to standard output */
    void printReport();

private:
    const MPI_Comm comm_;
    const std::string name_;
    std::stack<std::string> agents_;
    std::map<std::string, Accumulator> agents_all_;
    size_t batch_count_;
};

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* PROFILER_H_LQ95Z7KM */
