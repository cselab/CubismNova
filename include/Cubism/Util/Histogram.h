// File       : Histogram.h
// Created    : Mon Dec 23 2019 10:57:51 AM (+0100)
// Author     : Fabian Wermelinger
// Description: Histogram profiling on given MPI communicator
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef HISTOGRAM_H_6SPETKCI
#define HISTOGRAM_H_6SPETKCI

#include "Cubism/Common.h"
#include "Cubism/Util/Sampler.h"
#include <mpi.h>
#include <string>

NAMESPACE_BEGIN(Cubism)
/**
 * @addtogroup Util
 * @{ */
/** @brief Namespace for input/output operations
 *
 * @rst
 * The members of this namespace are optional utilities that are not required to
 * implement a working application.  Its components form the content of
 * ``libCubismUtil.a``.  An application using these utilities must link to
 * ``-lCubismUtil``.
 * @endrst
 */
NAMESPACE_BEGIN(Util)

/**
 * @ingroup Util MPI
 * @brief MPI profiling using histograms
 *
 * Collects samples for a profiled quantity of interest on individual ranks.
 * Can be used to detect inhomogeneities among MPI ranks.*/
class Histogram : public Sampler
{
public:
    /**
     * @brief Main histogram constructor
     * @param comm MPI communicator used for the profiling
     * @param name Name of the histogram
     * @param active Activator switch
     *
     * @rst
     * Extends the :ref:`sampler` class with MPI consolidation during
     * destruction of the object.  The consolidation generates a binary file
     * that can be post-processed using the :ref:`histbin` tool.  The activator
     * switch can be used to disable sample collection for large scale runs.
     * The size of the binary files can become large depending on the number of
     * ranks involved and the number of different quantities that are being
     * sampled.
     * @endrst
     */
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
/**  @} */
NAMESPACE_END(Cubism)

#endif /* HISTOGRAM_H_6SPETKCI */
