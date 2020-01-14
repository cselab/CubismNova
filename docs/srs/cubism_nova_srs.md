# Software Requirements Specification

The structure of this Software Requirements Specification (SRS) follows the
the [IEEE SRS template](https://github.com/rick4470/IEEE-SRS-Tempate).

## For CubismNova (alternatively Cubism-Nova)

```text
File       : cubism_nova_srs.md
Created    : Sun Nov 24 2019 02:35:12 PM (+0100)
Author     : Fabian Wermelinger (fabianw@mavt.ethz.ch)
             Michalis Chatzimanolakis (michaich@ethz.ch)
Version    : 1.0
Copyright 2019 ETH Zurich. All Rights Reserved.
```

Table of Contents
=================
<!-- vim-markdown-toc GFM -->

* [Revision History](#revision-history)
* [1 Introduction](#1-introduction)
  * [1.1 Purpose](#11-purpose)
  * [1.2 Document Conventions](#12-document-conventions)
  * [1.3 Intended Audience and Reading Suggestions](#13-intended-audience-and-reading-suggestions)
  * [1.4 Product Scope](#14-product-scope)
  * [1.5 References](#15-references)
* [2 Overall Description](#2-overall-description)
  * [2.1 Product Perspective](#21-product-perspective)
  * [2.2 Product Functions](#22-product-functions)
  * [2.3 User Classes and Characteristics](#23-user-classes-and-characteristics)
  * [2.4 Operating Environment](#24-operating-environment)
  * [2.5 Design and Implementation Constraints](#25-design-and-implementation-constraints)
  * [2.6 User Documentation](#26-user-documentation)
  * [2.7 Assumptions and Dependencies](#27-assumptions-and-dependencies)
* [3 External Interface Requirements](#3-external-interface-requirements)
  * [3.1 User Interfaces](#31-user-interfaces)
  * [3.2 Hardware Interfaces](#32-hardware-interfaces)
  * [3.3 Software Interfaces](#33-software-interfaces)
  * [3.4 Communications Interfaces](#34-communications-interfaces)
* [4 System Features](#4-system-features)
  * [4.1 System Feature 1](#41-system-feature-1)
  * [4.2 System Feature 2 (and so on)](#42-system-feature-2-and-so-on)
* [5 Other Non-Functional Requirements](#5-other-non-functional-requirements)
  * [5.1 Performance Requirements](#51-performance-requirements)
  * [5.2 Safety Requirements](#52-safety-requirements)
  * [5.3 Security Requirements](#53-security-requirements)
  * [5.4 Software Quality Attributes](#54-software-quality-attributes)
  * [5.5 Business Rules](#55-business-rules)
* [6 Other Requirements](#6-other-requirements)
  * [6.1 License](#61-license)
* [Appendix A: Glossary](#appendix-a-glossary)
* [Appendix B: Analysis Models](#appendix-b-analysis-models)
* [Appendix C: To Be Determined List](#appendix-c-to-be-determined-list)

<!-- vim-markdown-toc -->

## Revision History

| Name | Date | Reason For Changes | Version |
|------|------|--------------------|---------|
|      |      |                    |         |

## 1 Introduction

### 1.1 Purpose

This SRS describes the CubismNova library for structured uniform resolution or
(adaptive) multi-resolution (MR, AMR) resolution problems.

### 1.2 Document Conventions

**TBD**: Describe any standards or typographical conventions that were followed
when writing this SRS, such as fonts or highlighting that have special
significance. For example, state whether priorities for higher-level
requirements are assumed to be inherited by detailed requirements, or whether
every requirement statement is to have its own priority.

### 1.3 Intended Audience and Reading Suggestions

The document is written for the scientific community in the field of
computational and physical sciences, but may prove useful in other fields or
industrial applications as well.

**TODO**:

* Describe what the rest of the SRS contains
* How ist the SRS organized
* Suggest a reading order for each reader type (start with overview section and
proceed through sections that are most pertinent to each reader type)

### 1.4 Product Scope

CubismNova provides the fundamental framework for point-wise and stencil
based operations with a rich set of tools for basic operations such as
finite-differences or various interpolation schemes suitable and optimized for
recent High Performance Computing (HPC) architectures as well as smaller
workstations.  Example applications that make use of the CubismNova library (or
its predecessor Cubism) are:

* Compressible multicomponent flow solvers
* Incompressible multicomponent flow solvers
* Incompressible flow solvers containing solid obstacles coupled to machine
learning algorithms
* Lossy and loss-less high performance data compression algorithms

CubismNova is a completely revised and refactored C++ template library based on
the experience gathered by its successful predecessor Cubism.  On top of the
uniform structured grid principles inherited from Cubism, CubismNova extends
its capabilities with additional support for (Adaptive) Multi-Resolution in the
form of AMR, Block-Structured grids as well as Multigrid algorithms for various
data storage layouts.  The main goal of CubismNova is to separate the user from
the low-level optimization considerations embedded in the software.
User-friendliness and high-performance are contradictory concerns.  CubismNova
is redesigned to provide a meaningful interface to the user (application
designer) and at the same time maintain a high efficiency of hardware
utilization.

### 1.5 References

```
@article{hejazialhosseini2012a,
 address       = {Los Alamitos, CA, USA},
 author        = {Babak Hejazialhosseini and Diego Rossinelli and Christian Conti and Petros Koumoutsakos},
 doi           = {http://doi.ieeecomputersociety.org/10.1109/SC.2012.66},
 issn          = {2167-4329},
 journal       = {SC Conference},
 pages         = {1-12},
 publisher     = {IEEE Computer Society},
 title         = {High throughput software for direct numerical simulations of compressible two-phase flows},
 volume        = {0},
 year          = {2012}
}

@inproceedings{rossinelli2013a,
 acmid         = {2504565},
 address       = {New York, NY, USA},
 articleno     = {3},
 author        = {Rossinelli, Diego and Hejazialhosseini, Babak and Hadjidoukas, Panagiotis and Bekas, Costas and Curioni, Alessandro and Bertsch, Adam and Futral, Scott and Schmidt, Steffen J. and Adams, Nikolaus A. and Koumoutsakos, Petros},
 booktitle     = {Proceedings of the International Conference on High Performance Computing, Networking, Storage and Analysis},
 doi           = {10.1145/2503210.2504565},
 isbn          = {978-1-4503-2378-9},
 location      = {Denver, Colorado},
 numpages      = {13},
 pages         = {3:1--3:13},
 publisher     = {ACM},
 series        = {SC '13},
 title         = {11 {PFLOP/s} Simulations of Cloud Cavitation Collapse},
 url           = {http://doi.acm.org/10.1145/2503210.2504565},
 year          = {2013}
}

@inproceedings{hadjidoukas2015a,
 acmid         = {2820085},
 address       = {Edinburgh, Scotland, UK},
 author        = {Hadjidoukas, Panagiotis E. and Rossinelli, Diego and Hejazialhosseini, Babak and Koumoutsakos, Petros},
 booktitle     = {Proceedings of the 3rd International Conference on Exascale Applications and Software},
 isbn          = {978-0-9926615-1-9},
 location      = {Edinburgh, UK},
 numpages      = {6},
 pages         = {7--12},
 publisher     = {University of Edinburgh},
 series        = {EASC '15},
 title         = {From 11 to 14.4 {PFLOPs}: Performance Optimization for Finite Volume Flow Solver},
 url           = {http://dl.acm.org/citation.cfm?id=2820083.2820085},
 year          = {2015}
}

@inproceedings{hadjidoukas2015b,
 author        = {Panagiotis E. Hadjidoukas and Diego Rossinelli and Fabian Wermelinger and Jonas Sukys and Ursula Rasthofer and Christian Conti and Babak Hejazialhosseini and Petros Koumoutsakos},
 booktitle     = {Parallel Computing: On the Road to Exascale, Proceedings of the International Conference on Parallel Computing, ParCo 2015, 1-4 September 2015, Edinburgh, Scotland, {UK}},
 doi           = {10.3233/978-1-61499-621-7-767},
 pages         = {767--776},
 title         = {High throughput simulations of two-phase flows on {Blue Gene/Q}},
 url           = {http://dx.doi.org/10.3233/978-1-61499-621-7-767},
 year          = {2015}
}

@inproceedings{wermelinger2016a,
 acmid         = {2929914},
 address       = {New York, NY, USA},
 articleno     = {8},
 author        = {Wermelinger, Fabian and Hejazialhosseini, Babak and Hadjidoukas, Panagiotis and Rossinelli, Diego and Koumoutsakos, Petros},
 booktitle     = {Proceedings of the Platform for Advanced Scientific Computing Conference},
 doi           = {10.1145/2929908.2929914},
 numpages      = {10},
 pages         = {8:1--8:10},
 publisher     = {{ACM} Press},
 series        = {PASC '16},
 title         = {An Efficient Compressible Multicomponent Flow Solver for Heterogeneous {CPU/GPU} Architectures},
 url           = {http://doi.acm.org/10.1145/2929908.2929914},
 year          = {2016}
}

@article{wermelinger2018a,
 author        = {F. Wermelinger and U. Rasthofer and P.E. Hadjidoukas and P. Koumoutsakos},
 doi           = {10.1016/j.jocs.2018.01.008},
 journal       = {Journal of Computational Science},
 pages         = {217--225},
 publisher     = {Elsevier {BV}},
 title         = {Petascale simulations of compressible flows with interfaces},
 url           = {https://doi.org/10.1016%2Fj.jocs.2018.01.008},
 volume        = {26},
 year          = {2018}
}

@article{sukys2018a,
 author        = {Jonas {\v{S}}ukys and Ursula Rasthofer and Fabian Wermelinger and Panagiotis Hadjidoukas and Petros Koumoutsakos},
 doi           = {10.1137/17m1129684},
 journal       = {{SIAM} Journal on Scientific Computing},
 number        = {5},
 pages         = {B1361--B1390},
 publisher     = {Society for Industrial {\&} Applied Mathematics ({SIAM})},
 title         = {Multilevel Control Variates for Uncertainty Quantification in Simulations of Cloud Cavitation},
 url           = {https://doi.org/10.1137%2F17m1129684},
 volume        = {40},
 year          = {2018}
}

@article{verma2018a,
 author        = {Siddhartha Verma and Guido Novati and Petros Koumoutsakos},
 doi           = {10.1073/pnas.1800923115},
 journal       = {Proceedings of the National Academy of Sciences},
 number        = {23},
 pages         = {5849--5854},
 publisher     = {Proceedings of the National Academy of Sciences},
 title         = {Efficient collective swimming by harnessing vortices through deep reinforcement learning},
 url           = {https://doi.org/10.1073%2Fpnas.1800923115},
 volume        = {115},
 year          = {2018}
}

@article{rasthofer2019a,
 author        = {Rasthofer, U. and Wermelinger, F. and Karnakov, P. and {\v{S}}ukys, J. and Koumoutsakos, P.},
 doi           = {10.1103/PhysRevFluids.4.063602},
 issue         = {6},
 journal       = {Phys. Rev. Fluids},
 numpages      = {30},
 pages         = {063602},
 publisher     = {American Physical Society},
 title         = {Computational study of the collapse of a cloud with $12500$ gas bubbles in a liquid},
 url           = {https://link.aps.org/doi/10.1103/PhysRevFluids.4.063602},
 volume        = {4},
 year          = {2019}
}

@inproceedings{karnakov2019a,
 author        = {Petr Karnakov and Fabian Wermelinger and Michail Chatzimanolakis and Sergey Litvinov and Petros Koumoutsakos},
 booktitle     = {Proceedings of the Platform for Advanced Scientific Computing Conference on - {PASC} {\textquotesingle}19},
 doi           = {10.1145/3324989.3325727},
 publisher     = {{ACM} Press},
 title         = {A High Performance Computing Framework for Multiphase, Turbulent Flows on Structured Grids},
 url           = {https://doi.org/10.1145%2F3324989.3325727},
 year          = {2019}
}
```

## 2 Overall Description

### 2.1 Product Perspective
### 2.2 Product Functions
### 2.3 User Classes and Characteristics
### 2.4 Operating Environment
### 2.5 Design and Implementation Constraints
### 2.6 User Documentation
### 2.7 Assumptions and Dependencies

## 3 External Interface Requirements

### 3.1 User Interfaces
### 3.2 Hardware Interfaces
### 3.3 Software Interfaces
### 3.4 Communications Interfaces

## 4 System Features

### 4.1 System Feature 1
### 4.2 System Feature 2 (and so on)

## 5 Other Non-Functional Requirements

### 5.1 Performance Requirements
### 5.2 Safety Requirements
### 5.3 Security Requirements
### 5.4 Software Quality Attributes
### 5.5 Business Rules

## 6 Other Requirements

### 6.1 License

CubismNova is licensed under the BSD 2-Clause simplified license.

```text
CubismNova -- HPC library for structured uniform and adaptive multi resolution
Copyright © 2019 ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of ETH Zurich.
```

## Appendix A: Glossary
## Appendix B: Analysis Models
## Appendix C: To Be Determined List
