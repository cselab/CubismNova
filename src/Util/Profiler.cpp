// File       : Histogram.cpp
// Created    : Mon Dec 23 2019 11:02:20 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Histogram implementation
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Util/Profiler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

struct Stats {
    double nsamples; // use double for simplicity
    double mean;     // arithmetic mean
    double sdev;     // standard deviation
    double accu;     // accumulated value
    double vmin;     // minimum value
    double vmax;     // maximum value
    size_t rmin;     // index of minimum
    size_t rmax;     // index of maximum
};

Stats computeStats(const std::vector<double> &samples)
{
    Stats stat = {};

    stat.nsamples = samples.size();
    if (samples.size() < 2) {
        stat.mean = samples[0];
        stat.accu = stat.mean;
        stat.vmin = stat.mean;
        stat.vmax = stat.mean;
        return stat;
    }

    double sum = 0.0;
    double vmin = samples[0];
    double vmax = samples[0];
    size_t rmin = 0;
    size_t rmax = 0;
    for (size_t i = 0; i < samples.size(); ++i) {
        const double s = samples[i];
        sum += s;
        if (s < vmin) {
            vmin = s;
            rmin = i;
        }
        if (s > vmax) {
            vmax = s;
            rmax = i;
        }
    }
    const double mean = sum / samples.size();

    double sdev = 0.0;
    for (auto s : samples) {
        const double zmean = s - mean;
        sdev += zmean * zmean;
    }
    sdev = std::sqrt(sdev / (samples.size() - 1));

    stat.mean = mean;
    stat.sdev = sdev;
    stat.accu = sum;
    stat.vmin = vmin;
    stat.vmax = vmax;
    stat.rmin = rmin;
    stat.rmax = rmax;
    return stat;
}

void Profiler::printReport()
{
    int rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    double accumulated_time = 0.0;
    std::map<std::string, Stats> rank_stats;
    for (auto &agent : agents_all_) {
        const std::string &name = agent.first; // agent name
        auto &gs = agent.second;               // agent global rank stats

        const size_t id = std::hash<std::string>{}(name);
        size_t ref = id;
        MPI_Bcast(&ref, 1, MPI_UINT64_T, 0, comm_);
        if (id != ref) {
            throw std::runtime_error(
                "printReport: Non-symmetric agents across MPI ranks");
        }

        Stats s;
        const auto it = this->samples_.find(name);
        if (it != this->samples_.end()) {
            const std::vector<double> &data = it->second;

            // rank stats
            s = computeStats(data); // rank agent stats for this batch
            gs.batch_samples = data.size();
            gs.total_samples += gs.batch_samples;

            // send data to root
            std::vector<double> all_mean(size);
            MPI_Gather(&s.mean,
                       1,
                       MPI_DOUBLE,
                       all_mean.data(),
                       1,
                       MPI_DOUBLE,
                       0,
                       comm_);

            if (0 == rank) {
                const Stats rs = computeStats(all_mean); // rank statistics
                rank_stats[name] = rs;
                gs.total_time_mean += rs.mean;
                gs.total_time_accu += gs.batch_samples * rs.mean;
            } else {
                gs.total_time_mean += s.mean;
                gs.total_time_accu += s.accu;
            }
        }
        if (0 == rank) {
            accumulated_time += gs.total_time_accu;
        }
    }

    ++batch_count_;
    if (0 == rank) {
        std::string header = name_ + " profiler report: batch id = " +
                             std::to_string(batch_count_);
        std::transform(header.begin(), header.end(), header.begin(), ::toupper);
        std::printf("%s\n", header.c_str());

        // legend
        printf("  %-24s   %-10s %-10s min:%-10s:%-4s max:%-10s:%-4s %-7s "
               "-- %-10s %-10s %-8s %6s\n",
               "Name",
               "batch mean",
               "batch sdev",
               "value",
               "rank",
               "value",
               "rank",
               "samples",
               "mean",
               "total",
               "samples",
               "frac");
        for (const auto &agent : agents_all_) {
            const std::string &name = agent.first;
            const auto &gs = agent.second;

            Stats s;
            const auto it = this->samples_.find(name);
            if (it != this->samples_.end()) {
                s = rank_stats.find(name)->second;
            }
            printf(" [%-24s]: %.4e %.4e min:%.4e:%-4lu max:%.4e:%-4lu %-7lu "
                   "-- %.4e %.4e %-8lu %5.1f%%\n",
                   name.c_str(),
                   s.mean,
                   s.sdev,
                   s.vmin,
                   s.rmin,
                   s.vmax,
                   s.rmax,
                   gs.batch_samples,
                   gs.total_time_accu / gs.total_samples,
                   gs.total_time_accu,
                   gs.total_samples,
                   gs.total_time_accu / accumulated_time * 100.0);
        }
    }

    // clear this batch
    Sampler::clear();
}

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)
