// File       : Sampler.h
// Created    : Mon Dec 23 2019 10:54:43 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Sampler class used for Histogram and Profiler
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef SAMPLER_H_PV1MECI3
#define SAMPLER_H_PV1MECI3

#include "Cubism/Util/Timer.h"
#include <map>
#include <stack>
#include <string>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

/**
 * @ingroup Util
 * @brief Sample collector
 *
 * @rst
 * Used to collect time samples (by default) for a code section that is enclosed
 * by the ``seedSample()`` and ``collectSample()`` methods.  Used for profiling.
 * @endrst
 * */
class Sampler
{
public:
    /**
     * @brief Default constructor
     * @param active Activator switch
     */
    explicit Sampler(const bool active = true) : active_(active) {}
    virtual ~Sampler() {}

    using SampleMap = std::map<std::string, std::vector<double>>;

    /** @brief Seed new sample
     *
     * Pushes a new timer on the stack. */
    void seedSample()
    {
        if (active_) {
            timer_stack_.push(Timer());
        }
    }

    /**
     * @brief Collect a sample
     * @param name Name of the collected sample
     *
     * @rst
     * Pops the top timer off the stack and collects the measured time sample
     * for the given ``name``.
     * @endrst
     */
    void collectSample(const std::string &name)
    {
        if (active_ && !timer_stack_.empty()) {
            samples_[name].push_back(timer_stack_.top().stop());
            timer_stack_.pop();
        }
    }

    /**
     * @brief Pop a sample
     * @param name Name of the sample
     *
     * @rst
     * Pops the last (most recent) sample of ``name`` from the sample container.
     * @endrst
     */
    void popLast(const std::string &name)
    {
        if (active_ && !samples_[name].empty()) {
            samples_[name].pop_back();
        }
    }

    /**
     * @brief Append samples
     * @param s Sampler to take samples from
     *
     * @rst
     * Inserts all samples from ``s`` into this sampler.
     * @endrst
     */
    void append(const Sampler &s)
    {
        if (active_) {
            for (const auto &x : s.samples_) {
                auto &my_s = samples_[x.first];
                const auto &their_s = x.second;
                my_s.insert(my_s.end(), their_s.begin(), their_s.end());
            }
        }
    }

    /**
     * @brief Append single sample
     * @param name Name of the sample
     * @param samp Sample value
     *
     * @rst
     * Appends the sample ``samp`` to the list of samples with ``name``.
     * @endrst
     */
    void appendSample(const std::string &name, const double samp)
    {
        if (active_) {
            samples_[name].push_back(samp);
        }
    }

    /**
     * @brief Insert vector of samples
     * @param name Name of the samples
     * @param data Vector of samples to be inserted
     */
    void insert(const std::string &name, const std::vector<double> &data)
    {
        if (active_) {
            auto it = samples_.find(name);
            if (it == samples_.end()) {
                samples_[name] = data;
            } else {
                auto &my_s = it->second;
                my_s.insert(my_s.end(), data.begin(), data.end());
            }
        }
    }

    /**
     * @brief Perform addition with sample values
     * @param addto Name of samples to be added to
     * @param yours Vector of samples to be arithmetically added
     *
     * @rst
     * The number of samples in ``yours`` must be the same as the number of
     * samples in ``addto``, otherwise a runtime error is thrown.
     * @endrst
     */
    void addTo(const std::string &addto, const std::vector<double> &yours)
    {
        if (active_) {
            auto &mine = samples_.at(addto);
            if (mine.size() != yours.size()) {
                throw std::runtime_error("Add: vectors are of unequal length");
            }
            for (size_t i = 0; i < mine.size(); ++i) {
                mine[i] += yours[i];
            }
        }
    }

    /**
     * @brief Perform subtraction with sample values
     * @param from Name of samples to be subtracted from (minuend)
     * @param yours Vector of samples to be arithmetically subtracted
     * (subtrahend)
     *
     * @rst
     * The number of samples in ``yours`` must be the same as the number of
     * samples in ``from``, otherwise a runtime error is thrown.
     * @endrst
     */
    void subtractFrom(const std::string &from, const std::vector<double> &yours)
    {
        if (active_) {
            auto &mine = samples_.at(from);
            if (mine.size() != yours.size()) {
                throw std::runtime_error("Add: vectors are of unequal length");
            }
            for (size_t i = 0; i < mine.size(); ++i) {
                mine[i] -= yours[i];
            }
        }
    }

    /**
     * @brief Get the sample container
     * @return ``const`` reference to sample map
     */
    const SampleMap &getSamples() const { return samples_; }

    /** @brief Clear the sampler data
     *
     * After this method is called the sampler will contain zero samples. */
    void clear()
    {
        samples_.clear();
        std::stack<Timer>().swap(timer_stack_);
    }

protected:
    const bool active_;
    SampleMap samples_;
    std::stack<Timer> timer_stack_;
};

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* SAMPLER_H_PV1MECI3 */
