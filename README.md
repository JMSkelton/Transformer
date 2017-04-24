Transformer
===========

A Python library for transforming structures by performing atomic substitutions.

The package exposes a set of "convenience functions" for performing common operations such as generating solid solutions and antisite defects and batch-exporting results.
The "primitives" used to build these routines are factored into individual building blocks to allow for more advanced tasks be programmed quickly and easily.

The core algorithms make use of crystal symmetry (via `spglib`[1]) to reduce the computational workload, along with an optimised data layout (via `NumPy`[2]) for performance.

References
==========

[1] [https://atztogo.github.io/spglib/](https://atztogo.github.io/spglib/)
[2] [http://www.numpy.org/](http://www.numpy.org/)
