// File       : Timer.h
// Created    : Mon Dec 23 2019 11:07:58 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Simple timer
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef TIMER_H_1TOI9LUC
#define TIMER_H_1TOI9LUC

#include "Core/Common.h"

#include <chrono>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

class Timer
{
    using Clock = std::chrono::steady_clock;

public:
    Timer() : start_(clock_.now()) {}

    double getSeconds() const
    {
        return std::chrono::duration<double>(clock_.now() - start_).count();
    }

private:
    Clock clock_;
    Clock::time_point start_;
};

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* TIMER_H_1TOI9LUC */
