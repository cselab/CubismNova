.. File       : Field.rst
.. Created    : Thu Jan 16 2020 06:31:49 PM (+0100)
.. Author     : Fabian Wermelinger
.. Description: Block/Field.h documentation
.. Copyright 2020 ETH Zurich. All Rights Reserved.

.. _example-field:

Field
-----

Construct a block field
^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: Field_01.cpp
   :linenos:
   :language: cpp

Create a view to a block field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *view* is a cheap way to create a copy of another block field.  Copies for
views are *shallow* and therefore a view never owns memory.  A field view type
inherits from its underlying field type such that all of the field interface is
available for field views as well with different behavior for copy assignment.
A field view always performs shallow copy assignment.  This may not be desirable
if the user actually wants to copy data from a field to the field being viewed
at.  To enforce a deep copy the ``copyData()`` method can be used which must be
provided by the base field type.  A soft view provides a ``copy()`` method to
return a deep copy of the viewed field.

.. literalinclude:: Field_02.cpp
   :linenos:
   :language: cpp
