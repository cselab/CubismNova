// File       : Timer.h
// Created    : Mon Apr 01 2019 04:46:27 PM (+0200)
// Author     : Diego Rossinelli
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef TIMER_H_1TOI9LUC
#define TIMER_H_1TOI9LUC

#include "Common.h"
#include <sys/time.h>

NAMESPACE_BEGIN(Cubism)

class Timer
{
    struct timeval t_start, t_end;
    struct timezone t_zone;

public:
    void start() { gettimeofday(&t_start, &t_zone); }

    double stop()
    {
        gettimeofday(&t_end, &t_zone);
        return (t_end.tv_usec - t_start.tv_usec) * 1e-6 +
               (t_end.tv_sec - t_start.tv_sec);
    }
};

NAMESPACE_END(Cubism)

#endif /* TIMER_H_1TOI9LUC */
