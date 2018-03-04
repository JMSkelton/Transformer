# Transformer/Utilities/DevelopmentTools.py


# -------
# Imports
# -------

import numpy as np;

from Transformer import StructureSet;

from Transformer.Structure import Structure;

from Transformer.Utilities import MultiprocessingHelper;


# ---------
# Functions
# ---------

def MapStructures(parentStructure, structures1, structures2, tolerance = None, progressBar = True, useMP = False, mpNumProcesses = None):
    # If tolerance is not set, set it to the default from the Structure class.

    if tolerance == None:
        tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

    # Take the symmetry operations from the parent structure.

    parentSymmetryOperations = parentStructure.GetSymmetryOperations();

    # Map structures in the first result set to those in the second.

    inputItems = [
        (structure, structures2, parentSymmetryOperations, tolerance)
            for structure in structures1
        ];

    # The QueueMapFunction() routine in the MultiprocessingHelper module is used to parallelise over comaprisons of the reference structures in the first structure set to structures in the second.
    # If useMP is not set, QueueMapFunction() automatically falls back to a serial routine.

    structureMappings = MultiprocessingHelper.QueueMapFunction(
        _MapStructures_MapFunction, inputItems, maxNumProcesses = mpNumProcesses if useMP else 1, progressBar = progressBar
        );

    return structureMappings;

def _MapStructures_MapFunction(refStructure, compareStructures, parentSymmetryOperations, tolerance):
    refAtomCount = refStructure.GetAtomCount();

    # Range of atom indices to compare positions over.

    indexRanges = [(0, refAtomCount)];

    # Expand the reference structure with the parent symmetry operations.

    refTransformedPositions = None;

    # Cache symmetry expansions of the structures in compareStructures for efficiency.

    symmetryExpansionsCache = [
        None for _ in range(0, len(compareStructures))
        ];

    # Compare the transformed structures to each structure in the second result set.

    indices = [];

    for j, compareStructure in enumerate(compareStructures):
        # First, check the two structures have the same number of atoms.

        equivalent = compareStructure.GetAtomCount() == refAtomCount;

        # Depending on whether compareLatticeVectors/compareAtomTypes are set, compare the lattice vectors and/or atom-type numbers.

        if equivalent:
            latticeVectors1 = refStructure.GetLatticeVectorsNumPy(copy = False);
            latticeVectors2 = compareStructure.GetLatticeVectorsNumPy(copy = False);

            equivalent = np.all(np.abs(latticeVectors1 - latticeVectors2) < tolerance);

        if equivalent:
            atomTypeNumbers1 = refStructure.GetAtomTypeNumbersNumPy(copy = False);
            atomTypeNumbers2 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);

            equivalent = (atomTypeNumbers1 == atomTypeNumbers2).all();

        # Compare the atom positions.

        if equivalent:
            # Lazy initialisation of refTransformedPositions.

            if refTransformedPositions is None:
                refTransformedPositions = StructureSet._GenerateSymmetryExpandedPositions(
                    refStructure, parentSymmetryOperations, tolerance
                    );

            # Compare the atom positions in compareStructure to the set of symmetry-expanded reference positions.

            equivalent = StructureSet._CompareAtomPositions(
                compareStructure.GetAtomPositionsNumPy(copy = False), refTransformedPositions, tolerance, indexRanges
                );

            if not equivalent:
                # If we have already generated the symmetry-expanded positions for compareStructure, take them from the cache.
                # If not, generate and cache them.

                if symmetryExpansionsCache[j] is None:
                    symmetryExpansionsCache[j] = StructureSet._GenerateSymmetryExpandedPositions(
                        compareStructure, parentSymmetryOperations, tolerance
                        );

                compareTransformedPositions = symmetryExpansionsCache[j];

                # Compare each set of symmetry-expanded positions for compareStructure with the symmetry-expanded reference positions.

                for comparePositions in compareTransformedPositions:
                    equivalent = StructureSet._CompareAtomPositions(
                        comparePositions, refTransformedPositions, tolerance, indexRanges
                        );

                    if equivalent:
                        break;

        if equivalent:
            indices.append(j);

    # Return the list of indices.

    return indices;
