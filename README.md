Transformer
===========

`Transformer` is a Python library for transforming structures by performing atomic substitutions.

The package exposes a framework of "convenience functions" for performing common operations, including generating solid solutions and antisite defects, and batch-importing/exporting result sets.

The logic used to build up these routines is separated into lower-level "primitives" that can be used to program more advanced tasks.

Installation
============

The `Transformer` library works in Python >= 2.7 and Python 3 (so far tested on macOS and Debian Linux 8.8), and just needs to be placed in the package search path to work (e.g. add to `PYTHONPATH`).

The library requires the `NumPy` ([Ref. 1](#Ref1)) and `spglib` ([Ref. 2](#Ref2)) packages; both can be installed using `pip` on Linux/macOS, or built from source (see links to code documentation).

Installing the `Cython` ([Ref. 3](#Ref3)) and `tqdm` ([Ref. 4](#Ref4)) packages (both also available *via* `pip`) enables optional performance enhancements and "eye candy".
`Cython` in particular is highly recommended, as it allows core routines to be replaced by optimised C code and leads to substantial performance enhancements.

Examples
========

At present, `Transformer` is not fully documented.
The best way to learn how to use it is to look at the [Examples](./Examples) folder, which contains a few examples illustrating how to perform common functions:

1. [SnS/Se solid solutions](./Examples/Example_SnS-Se-SolidSolution.py): Evaluate the full set of inequivalent structures in an S/Se solid solution built from a 2&times;1&times;2 supercell of *Pnma* SnS with 32 atoms.

2. [Antisite defects in Cu<sub>2</sub>ZnSnS<sub>4</sub>](./Examples/Example_CZTS-AntisiteDefects.py): Generate configurations with 1-2 antisite defects in a 2&times;2&times;1 supercell of kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub>.

3. [ZnS "cation mutation"](./Examples/Example_ZnS-CationMutation.py): Starting from a 1&times;1&times;2 supercell expansion of the conventional ZnS cell, generate CuGaS<sub>2</sub> and Stannite/Kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub> by successive cation substitutions.

4. [Cu<sub>2</sub>ZnSnS<sub>4</sub> vacancies](./Examples/Example_CZTS-Vacancies-AIMS.py): Generate 2&times;2&times;1 supercells of kesterite Cu<sub>2</sub>ZnSnS<sub>4</sub> with Cu, Zn, Sn and S defects, and cation vacancies with compensating anion vacancies (i.e. Schottky defects).
   This example imports and exports FHI-AIMS rather than VASP POSCAR files.

Working TODO List
=================

* Finish adding docstrings (so far only the `Structure` and `StructureSet` modules have been documented).
* Add examples for using the multithreaded substitution routine and `Filters` framework.
* Add examples for using the `Screening` module.
* Implement statistical-thermodynamics analysis of data from structure sets generated by `Transformer`.
* Implement functionality for identifying and working with interstitials.

References
==========

1. <a name="Ref1"></a>[https://atztogo.github.io/spglib/](https://atztogo.github.io/spglib/)
2. <a name="Ref2"></a>[http://www.numpy.org/](http://www.numpy.org/)
3. <a name="Ref3"></a>[http://cython.org/](http://cython.org/)
4. <a name="Ref4"></a>[https://pypi.python.org/pypi/tqdm](https://pypi.python.org/pypi/tqdm); [https://github.com/noamraph/tqdm](https://github.com/noamraph/tqdm)
