# Transformer/Utilities/DevelopmentTools.py by J. M. Skelton


# -------
# Imports
# -------

import numpy as np;

from Transformer import StructureSet;

# If available, import the Cython-optimised merging routines from the _StructureSet module.

_Cython = False;

try:
    from Transformer import _StructureSet;

    _Cython = True;
except ImportError:
    pass;

from Transformer.Structure import Structure;


# ---------
# Functions
# ---------

def MapResultSetStructures(parentStructure, structures1, structures2, compareLatticeVectors = True, compareAtomTypes = True, tolerance = None):
    # If tolerance is not set, set it to the default from the Structure class.

    if tolerance == None:
        tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

    # Take the symmetry operations from the parent structure.

    parentSymmetryOperations = parentStructure.GetSymmetryOperations();

    # Map structures in the first result set to those in the second.

    structureMappings = [];

    for i, refStructure in enumerate(structures1):
        refAtomCount = refStructure.GetAtomCount();

        # Range of atom indices to compare positions over.

        indexRange = [(0, refAtomCount)];

        # Expand the reference structure to a set of comparison structures by applying the parent symmetry operations.

        comparePositions = None;

        # Compare the transformed structures to each structure in the second result set.

        indices = [];

        for j, compareStructure in enumerate(structures2):
            # First, check the two structures have the same number of atoms.

            equivalent = compareStructure.GetAtomCount() == refAtomCount;

            # Depending on whether compareLatticeVectors/compareAtomTypes are set, compare the lattice vectors and/or atom-type numbers.

            if equivalent and compareLatticeVectors:
                equivalent = refStructure.CompareLatticeVectors(compareStructure, tolerance = tolerance);

            if equivalent and compareAtomTypes:
                equivalent = refStructure.CompareAtomTypeNumbers(compareStructure);

            # Compare the atom positions in structure to the set of symmetry-transformed positions in comparePositions.

            if equivalent:
                # Lazy initialisation of comparePositions.

                if comparePositions is None:
                    # If possible, use the Cython-optimised symmetry-transformation routine.

                    if _Cython:
                        comparePositions = _StructureSet._GenerateSymmetryTransformedPositions(
                            refStructure, parentSymmetryOperations, tolerance
                            );
                    else:
                        comparePositions = StructureSet._GenerateSymmetryTransformedPositions(
                            refStructure, parentSymmetryOperations, tolerance
                            );

                if _Cython:
                    # If possible, use the Cython-optimised comparison routine.

                    equivalent = _StructureSet._CompareAtomPositions(
                        compareStructure.GetAtomPositionsNumPy(copy = False), comparePositions, tolerance, indexRange
                        );
                else:
                    equivalent = StructureSet._CompareAtomPositions(
                        compareStructure.GetAtomPositionsNumPy(copy = False), comparePositions, tolerance, indexRange
                        );

            if equivalent:
                indices.append(j);

        structureMappings.append(indices);

    return structureMappings;
