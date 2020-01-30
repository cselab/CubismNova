.. File       : CartesianHDF.rst
.. Created    : Thu Jan 30 2020 06:04:19 PM (+0100)
.. Author     : Fabian Wermelinger
.. Description: IO/CartesianHDF.h documentation
.. Copyright 2020 ETH Zurich. All Rights Reserved.

CartesianHDF.h
----------

.. doxygenfunction:: CartesianWriteHDF(const std::string&, const std::string&, const Grid&, const Mesh&, const double, const Dir, const bool)
   :project: CubismNova

.. doxygenfunction:: CartesianWriteHDF(const std::string&, const std::string&, const Grid&, const double, const Dir, const bool)
   :project: CubismNova

.. doxygenfunction:: CartesianReadHDF(const std::string&, Grid&, const Mesh&, const Dir)
   :project: CubismNova

.. doxygenfunction:: CartesianReadHDF(const std::string&, Grid&, const Dir)
   :project: CubismNova
