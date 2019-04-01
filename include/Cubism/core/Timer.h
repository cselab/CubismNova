// File       : Timer.h
// Created    : Mon Apr 01 2019 04:46:27 PM (+0200)
// Author     : Diego Rossinelli
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef TIMER_H_1TOI9LUC
#define TIMER_H_1TOI9LUC

#include <sys/time.h>

#include "Common.h"

NAMESPACE_BEGIN(cubism)

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

NAMESPACE_END(cubism)

#endif /* TIMER_H_1TOI9LUC */
