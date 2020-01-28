// File       : Histogram.cpp
// Created    : Mon Dec 23 2019 11:02:20 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Histogram implementation
// Copyright 2019 ETH Zurich. All Rights Reserved.

#include "Cubism/Util/Histogram.h"
#include <cmath>
#include <stdexcept>
#include <utility>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

struct FileHeader {
    int nranks;  // number of ranks
    int nstats;  // number of stats per rank
    int nfloats; // number of 64bit floating point numbers per stat
};

struct Stats {
    double nsamples; // use double for simplicity
    double mean;     // arithmetic mean
    double sdev;     // standard deviation
    double accu;     // accumulated value
    double vmin;     // minimum value
    double vmax;     // maximum value
};

static std::map<std::string, Stats>
computeStats(const typename Sampler::SampleMap &samples)
{
    std::map<std::string, Stats> stats;
    for (auto &key : samples) {
        Stats stat = {};
        const std::vector<double> &data = key.second;
        stat.nsamples = data.size();
        if (data.size() < 2) {
            stats[key.first] = stat;
            continue;
        }
        double sum = 0.0;
        double vmin = data[0];
        double vmax = data[0];
        for (auto s : data) {
            sum += s;
            vmin = (s < vmin) ? s : vmin;
            vmax = (s > vmax) ? s : vmax;
        }
        const double mean = sum / data.size();

        double sdev = 0.0;
        for (auto s : data) {
            const double zmean = s - mean;
            sdev += zmean * zmean;
        }
        sdev = std::sqrt(sdev / (data.size() - 1));

        stat.mean = mean;
        stat.sdev = sdev;
        stat.accu = sum;
        stat.vmin = vmin;
        stat.vmax = vmax;
        stats[key.first] = stat;
    }
    return stats;
}

void Histogram::consolidate_()
{
    int rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    homogenizeCollection_();

    std::vector<size_t> char_bytes;
    std::vector<std::string> keys;
    size_t sum_char = 0;
    for (auto &key : samples_) {
        const size_t s = key.first.length();
        sum_char += s;
        char_bytes.push_back(s);
        keys.push_back(key.first);
    }

    const size_t nvalues = sizeof(Stats) / sizeof(double);
    const auto stats = computeStats(samples_);

    FileHeader header;
    header.nranks = size;
    header.nstats = static_cast<int>(stats.size());
    header.nfloats = nvalues;

    const size_t data_start =
        sizeof(FileHeader) + char_bytes.size() * sizeof(int) + sum_char;
    const size_t my_offset =
        data_start + rank * stats.size() * nvalues * sizeof(double);

    MPI_File fh;
    MPI_Status st;
    std::string name = "hist_" + name_ + ".bin";
    if (MPI_File_open(comm_,
                      name.c_str(),
                      MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN,
                      MPI_INFO_NULL,
                      &fh)) {
        throw std::runtime_error("Consolidate: Can not open MPI file");
    }

    if (0 == rank) {
        // header and keys
        MPI_File_seek(fh, 0, MPI_SEEK_SET);
        MPI_File_write(fh, &header, sizeof(FileHeader), MPI_BYTE, &st);
        for (size_t i = 0; i < samples_.size(); ++i) {
            int cb = static_cast<int>(char_bytes[i]);
            MPI_File_write(fh, &cb, 1, MPI_INT, &st);
            MPI_File_write(fh, keys[i].c_str(), cb, MPI_CHAR, &st);
        }
    } else {
        MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
    }
    for (const auto &s : stats) {
        // data
        MPI_File_write(fh, &s.second, sizeof(Stats), MPI_BYTE, &st);
    }
    MPI_File_close(&fh);
}

void Histogram::homogenizeCollection_()
{
    // if sample collection differs among participating ranks, homogenize
    // missing samples in current collection
    int rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    std::vector<char> cstream;
    int nchar = 0;
    for (auto &k : samples_) {
        for (const char *c = k.first.c_str(); *c; ++c, ++nchar) {
            cstream.push_back(*c);
        }
        cstream.push_back('\0');
        ++nchar;
    }

    std::vector<int> all_sizes(size);
    MPI_Allgather(&nchar, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, comm_);

    int sum_sizes = 0;
    std::vector<int> all_offsets = all_sizes;
    const int s0 = all_offsets[0];
    for (size_t i = 0; i < all_sizes.size(); ++i) {
        sum_sizes += all_sizes[i];
        all_offsets[i] -= s0;
    }
    std::vector<char> all_char(sum_sizes);
    MPI_Allgatherv(cstream.data(),
                   nchar,
                   MPI_CHAR,
                   all_char.data(),
                   all_sizes.data(),
                   all_offsets.data(),
                   MPI_CHAR,
                   comm_);

    const char *c = all_char.data();
    const std::vector<double> empty(0);
    while (sum_sizes > 0) {
        const std::string n(c);
        if (n != "") {
            insert(n, empty);
        }
        c += (n.length() + 1);
        sum_sizes -= (n.length() + 1);
    }
}

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)
