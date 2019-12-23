// File       : Sampler.h
// Created    : Mon Dec 23 2019 10:54:43 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Sampler class used for Histogram and Profiler
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef SAMPLER_H_PV1MECI3
#define SAMPLER_H_PV1MECI3

#include "Timer.h"

#include <map>
#include <stack>
#include <string>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

class Sampler
{
public:
    Sampler(const bool active = true) : active_(active) {}
    virtual ~Sampler() {}

    using SampleMap = std::map<std::string, std::vector<double>>;

    void seedSample()
    {
        if (active_) {
            timer_stack_.push(Timer());
        }
    }

    void collectSample(const std::string &name)
    {
        if (active_ && !timer_stack_.empty()) {
            samples_[name].push_back(timer_stack_.top().getSeconds());
            timer_stack_.pop();
        }
    }

    void popLast(const std::string &name)
    {
        if (active_ && !samples_[name].empty()) {
            samples_[name].pop_back();
        }
    }

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

    void appendSample(const std::string &name, const double samp)
    {
        if (active_) {
            samples_[name].push_back(samp);
        }
    }

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

    const SampleMap &getSamples() const { return samples_; }

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
