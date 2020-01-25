.. File       : Field.rst
.. Created    : Thu Jan 16 2020 06:31:49 PM (+0100)
.. Author     : Fabian Wermelinger
.. Description: Block/Field.h documentation
.. Copyright 2020 ETH Zurich. All Rights Reserved.

.. _field:

Field
-----

This header defines various block field types.  The essential block field is a
*scalar* field.  Components of higher rank tensors are built using the
corresponding scalar field type.  The memory layout of a field is a structure of
arrays.  Face fields are handled separately for each dimension in ``DIM``.  The
``FaceContainer`` type is a convenience type for a compound of ``DIM`` face
fields with the corresponding index range for each direction.

.. doxygenstruct:: FieldState
   :project: CubismNova
   :members:

.. doxygenclass:: Field
   :project: CubismNova
   :members:

.. doxygenclass:: TensorField
   :project: CubismNova
   :members:

.. doxygenclass:: FaceContainer
   :project: CubismNova
   :members:

.. doxygenclass:: FieldContainer
   :project: CubismNova
   :members:

.. doxygenclass:: FieldView
   :project: CubismNova
   :members:

.. doxygentypedef:: CellField
   :project: CubismNova

.. doxygentypedef:: NodeField
   :project: CubismNova

.. doxygentypedef:: FaceField
   :project: CubismNova

.. doxygentypedef:: VectorField
   :project: CubismNova
