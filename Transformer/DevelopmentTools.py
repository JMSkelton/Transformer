# Transformer/DevelopmentTools.py by J. M. Skelton


# -------
# Imports
# -------

import numpy as np;

from Transformer import StructureTools;

# If available, import the Cython-optimised merging routines from the _StructureTools module.

if StructureTools.CythonImports:
    from Transformer import _StructureTools;

from Transformer import Structure;


# ---------
# Functions
# ---------

def MapResultSetStructures(parentStructure, structures1, structures2, compareLatticeVectors = True, compareAtomTypes = True, tolerance = None):
    # If tolerance is not set, set it to the default from the Structure class.

    if tolerance == None:
        tolerance = Structure.Structure.DefaultSymmetryEquivalenceTolerance;

    # Take the symmetry operations from the parent structure.

    parentSymmetryOperations = parentStructure.GetSymmetryOperations();

    # Map structures in the first result set to those in the second.

    structureMappings = [];

    for i, refStructure in enumerate(structures1):
        # Expand the reference structure to a set of comparison structures by applying the parent symmetry operations.

        comparePositions = None;

        # Compare the transformed structures to each structure in the second result set.

        indices = [];

        for j, compareStructure in enumerate(structures2):
            equivalent = True;

            # Depending on whether compareLatticeVectors/compareAtomTypes are set, compare the lattice vectors and/or atom-type numbers.

            if compareLatticeVectors:
                equivalent = refStructure.CompareLatticeVectors(compareStructure, tolerance = tolerance);

            if equivalent and compareAtomTypes:
                equivalent = refStructure.CompareAtomTypeNumbers(compareStructure);

            # Compare the atom positions in structure to the set of symmetry-transformed positions in comparePositions.

            if equivalent:
                # Lazy initialisation of comparePositions.

                if comparePositions is None:
                    # If possible, use the Cython-optimised symmetry-transformation routine.

                    if StructureTools.CythonImports:
                        comparePositions = _StructureTools._MergeStructureSet_GenerateSymmetryTransformedPositions(
                            refStructure, parentSymmetryOperations
                            );
                    else:
                        comparePositions = StructureTools._MergeStructureSet_GenerateSymmetryTransformedPositions(
                            refStructure, parentSymmetryOperations
                            );

                if StructureTools.CythonImports:
                    # If possible, use the Cython-optimised comparison routine.

                    equivalent = _StructureTools._MergeStructureSet_ComparePositions(
                        compareStructure.GetAtomPositionsNumPy(copy = False), comparePositions, tolerance
                        );
                else:
                    compareResult = np.all(
                        np.abs(comparePositions - compareStructure.GetAtomPositionsNumPy(copy = False)) < tolerance, axis = (1, 2)
                        );

                    equivalent = compareResult.any();

            if equivalent:
                indices.append(j);

        structureMappings.append(indices);

    return structureMappings;
