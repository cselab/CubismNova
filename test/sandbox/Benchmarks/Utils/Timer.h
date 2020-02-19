// File       : Timer.h
// Created    : Mon Dec 23 2019 11:07:58 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Simple timer
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef TIMER_H_QBP9K0C1
#define TIMER_H_QBP9K0C1

#include <chrono>

namespace Utils
{
class Timer
{
    using Clock = std::chrono::steady_clock;

public:
    Timer() : start_(clock_.now()) {}

    void start() { start_ = clock_.now(); }
    double stop() const
    {
        return std::chrono::duration<double>(clock_.now() - start_).count();
    }

private:
    Clock clock_;
    Clock::time_point start_;
};
} // namespace Utils

#endif /* TIMER_H_QBP9K0C1 */
