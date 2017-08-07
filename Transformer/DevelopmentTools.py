# Transformer/DevelopmentTools.py by J. M. Skelton


# -------
# Imports
# -------

import numpy as np;

from Transformer import StructureTools;
from Transformer import ConvenienceFunctions;

from Transformer.Structure import Structure;


# ---------
# Functions
# ---------

def MapResultSetStructures(parentStructure, structures1, structures2, compareLatticeVectors = True, compareAtomTypes = True, tolerance = None):
    # If tolerance is not set, set it to the default from the Structure class.

    if tolerance == None:
        tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

    # Take the symmetry operations from the parent structure.

    symmetryOperations = parentStructure.GetSymmetryOperations();

    # Map structures in the first result set to those in the second.

    structureMappings = [];

    for i, structureRef in enumerate(structures1):
        # Expand the reference structure to a set of comparison structures by applying the parent symmetry operations.

        transformedStructures = StructureTools.GenerateSymmetryTransformedStructures(
            structureRef, symmetryOperations
            );

        comparePositions = transformedStructures.view(dtype = np.float64).reshape((len(transformedStructures), structureRef.GetAtomCount(), 4))[:, :, 1:];

        # Compare the transformed structures to each structure in the second result set.

        indices = [];

        for j, structure in enumerate(structures2):
            equivalent = True;

            # Depending on whether compareLatticeVectors/compareAtomTypes are set, compare the lattice vectors and/or atom-type numbers.

            if compareLatticeVectors:
                equivalent = structureRef.CompareLatticeVectors(structure, tolerance = tolerance);

            if equivalent and compareAtomTypes:
                equivalent = structureRef.CompareAtomTypeNumbers(structure);

            # Compare the atom positions in structure to the set of symmetry-transformed positions in comparePositions.

            if equivalent:
                compareResult = np.all(
                    np.abs(comparePositions - structure.GetAtomPositionsNumPy(copy = False)) < tolerance, axis = (1, 2)
                    );

                equivalent = compareResult.any();

            if equivalent:
                indices.append(j);

        structureMappings.append(indices);

    return structureMappings;
