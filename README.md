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

1. [SnS/Se solid solutions](./Examples/Example_SnS-Se-SolidSolution.py): Evaluate the full set of inequivalent structures in an S/Se solid solution built from a 2&times;1&times;2 supercell of *Pnma* SnS with 32 atoms.

2. [Antisite defects in Cu<sub>2</sub>ZnSnS<sub>4</sub>](./Examples/Example_CZTS-AntisiteDefects.py): Generate configurations with 1-4 antisite defects in a 2&times;2&times;1 supercell of kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub>.

3. [ZnS "cation mutation"](./Examples/Example_ZnS-CationMutation.py): Starting from a 1&times;1&times;2 supercell expansion of the conventional ZnS cell, generate CuGaS<sub>2</sub> and Stannite/Kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub> by successive cation substitutions.

4. [Cu<sub>2</sub>ZnSnS<sub>4</sub> vacancies](./Examples/Example_CZTS-Vacancies-AIMS.py): Generate 2&times;2&times;1 supercells of kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub> with Cu, Zn, Sn and S defects, and cation vacancies with compensating anion vacancies (i.e. Schottky defects).
   This example uses the FHI-AIMS I/O routines.

There are also a set of more advanved "development tests" in the [Development Tests](./DevelopmentTests) folder:

1. [Solid solutions "shortcut" check](./DevelopmentTests/DevelopmentTest_SolidSolution-ShortcutCheck.py): Runs the SnS/Se solid-solutions example with and without the "shortcut", then reads back the result sets and verifies that they are equivalent.

TODO List
=========

* Add `multiprocessing` support.
* Add a convenience function for Schottky defects.
* Add functionality for identifying interstitial sites for Frenkel defects.

References
==========

1. <a name="Ref1"></a>[https://atztogo.github.io/spglib/](https://atztogo.github.io/spglib/)
2. <a name="Ref2"></a>[http://www.numpy.org/](http://www.numpy.org/)
