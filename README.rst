.. File       : README.rst
.. Created    : Tue Jan 14 2020 06:34:44 PM (+0100)
.. Author     : Fabian Wermelinger
.. Description: CubismNova main README file
.. Copyright 2020 ETH Zurich. All Rights Reserved.

**********
CubismNova
**********

.. image:: https://img.shields.io/badge/version-0.0.1-blue
   :alt: Version
.. image:: https://img.shields.io/badge/license-BSD-green
   :target: LICENSE
   :alt: License
.. image:: https://circleci.com/gh/cselab/CubismNova.svg?style=shield
   :target: https://circleci.com/gh/cselab/CubismNova
   :alt: Build Status

Introduction
************

CubismNova is a C++ template library used for solving Partial Differential
Equations (PDEs) on structured uniform or stretched grids as well as
block-structured or adaptively refined grids (AMR).  The library provides data
structures for point-wise and stencil operations with support for efficient halo
(ghost cell) communication using the Message Passing Interface (MPI).  A toolbox
with vectorized kernels for common operations such as finite difference
operators, WENO reconstruction, interpolation, restriction and prolongation
operators as well as various data compression schemes is available to the user.
Extended data structures for use with various time integration schemes is
available as well.

The library is a full refactoring of its successful predecessor Cubism that has
won the Gordon Bell award for a compressible multicomponent flow problem in 2013
:cite:`rossinelli2013a`.  Further optimizations on the same code are presented
in :cite:`hadjidoukas2015a` and :cite:`hadjidoukas2015b`.  The refactored
library offers easier access for the community by separating high performance
computing (HPC) concepts from the user of the library.  The library user must be
concerned with the algorithm design depending on the problem that needs to be
solved.  The refactored library further offers integrated multigrid solvers and
compression algorithms to reduce the I/O overhead at scale.  Moreover, the
refactored library takes into account suitable data structures for use with
heterogeneous accelerators :cite:`wermelinger2016a`.  Apart from compressible
multicomponent flow simulations (:cite:`sukys2018a`, :cite:`wermelinger2018a`,
:cite:`rasthofer2019a`), the library is also used for incompressible multi-phase
flow (:cite:`karnakov2019a`) as well as incompressible flow with collective
swimmers (:cite:`verma2018a`).

Installation
************

CubismNova can be downloaded from GitHub_:

.. code:: console

   $ git clone --recurse-submodules https://github.com/cselab/CubismNova
   $ cd CubismNova

The library can be compiled using ``cmake``.  A working MPI implementation is
required and the ``mpicc`` and ``mpic++`` compiler wrappers must be in the
``PATH`` environment variable.  A debug build can be generated with

.. code:: console

   $ ./cmake_init.sh debug <install path>
   $ cd debug
   $ make -j && make test
   $ make install
   $ cd .. && rm -rf debug

This assumes that your starting directory is the project root of ``CubismNova``.
An optimized build is likewise generated with

.. code:: console

   $ ./cmake_init.sh release <install path>
   $ cd release
   $ make -j && make test
   $ make install
   $ cd .. && rm -rf release

Instead of ``release`` you can use any other token except ``debug``.  If the
``<insall path>`` is a system directory use ``sudo make install`` instead.

License
*******

`BSD`_ © 2019 ETH Zurich

.. _BSD: LICENSE
.. _GitHub: https://github.com/cselab/CubismNova
.. bibliography:: bibtex/references.bib
