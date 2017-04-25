Transformer
===========

`Transformer` is a Python library for transforming structures by performing atomic substitutions.

The package exposes a set of "convenience functions" for performing common operations, including generating solid solutions and antisite defects, and batch-exporting result sets.

The logic used to build up these routines is separated into lower-level "primitives" that can be used to program more advanced tasks.

The routines make use of crystal symmetry (via `spglib` - [Ref. 1](#Ref1)) to avoid, at least in most cases, computing and storing large sets of intermediate structures.
The core `Structure` object uses an optimised data layout (via `NumPy` - [Ref. 2](#Ref2)) for efficiency.

Installation
============

The `Transformer` library works in Python >= 2.7 and Python 3 (tested on macOS), and just needs to be placed in the package search path to work (e.g. add to `PYTHONPATH`).

The library requires the `NumPy` and `spglib` packages; both can be installed using `pip` on Linux/macOS, or built from source (see links to code documentation in [Ref. 1](#Ref1) and [Ref. 2](#Ref2)).

Examples
========

At present, apart from inline comments in the code files, `Transformer` has not been fully documented.
However, the [Examples](./Examples) folder contains a few examples illustrating how to perform common functions:

1. [Antisite defects in Cu<sub>2</sub>ZnSnS<sub>4</sub>](./Examples/Example_CZTS-AntisiteDefects.py): Generate configurations with 1-4 antisite defects in a 2&times;2&times;1 supercell of kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub>.

2. [ZnS "cation mutation"](./Examples/Example_ZnS-CationMutation.py): Starting from a 1&times;1&times;2 supercell expansion of the conventional ZnS cell, generate CuGaS<sub>2</sub> and Stannite/Kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub> by successive cation substitutions.

References
==========

1. <a name="Ref1"></a>[https://atztogo.github.io/spglib/](https://atztogo.github.io/spglib/)
2. <a name="Ref2"></a>[http://www.numpy.org/](http://www.numpy.org/)
