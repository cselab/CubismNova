.. File       : README.rst
.. Created    : Tue Jan 14 2020 06:34:44 PM (+0100)
.. Author     : Fabian Wermelinger
.. Description: CubismNova main README file
.. Copyright 2020 ETH Zurich. All Rights Reserved.

**********
CubismNova
**********

.. image:: https://circleci.com/gh/cselab/CubismNova.svg?style=shield
   :target: https://circleci.com/gh/cselab/CubismNova
   :alt: Build Status
.. image:: https://readthedocs.org/projects/cubismnova/badge/?version=latest
   :target: https://cubismnova.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://codecov.io/gh/cselab/CubismNova/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/cselab/CubismNova
   :alt: Code Coverage

Introduction
============

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

The library is a full refactoring of its successful predecessor Cubism
:cite:t:`hejazialhosseini2012a` that has later won the Gordon Bell award for a
compressible multicomponent flow problem in 2013 :cite:t:`rossinelli2013a`.
Further optimizations on the same code are presented in
:cite:t:`hadjidoukas2015a` and :cite:t:`hadjidoukas2015b`.  The refactored
library offers easier access for the community by separating high performance
computing (HPC) concepts from the user of the library.  The library user must be
concerned with the algorithm design depending on the problem that needs to be
solved.  The refactored library further offers integrated multigrid solvers and
compression algorithms to reduce the I/O overhead at scale.  Moreover, the
refactored library takes into account suitable data structures for use with
heterogeneous accelerators :cite:t:`wermelinger2016a`.  Apart from compressible
multicomponent flow simulations (:cite:t:`sukys2018a`,
:cite:t:`wermelinger2018a`, :cite:t:`rasthofer2019a`), the library is also used
for incompressible multi-phase flow (:cite:t:`karnakov2019a`) as well as
incompressible flow with collective swimmers (:cite:t:`verma2018a`).

Documentation
=============

The software documentation is hosted at https://cubismnova.readthedocs.io

Installation
============

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

Versioning
==========

This software follows the `semantic versioning specification`_.

License
=======

`BSD`_ © 2019 ETH Zurich

.. _BSD: LICENSE
.. _GitHub: https://github.com/cselab/CubismNova
.. _semantic versioning specification: https://semver.org/

References
==========

.. bibliography::
