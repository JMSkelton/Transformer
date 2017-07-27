# Transformer/StructureTools.py by J. M. Skelton


# -------
# Imports
# -------

import warnings;

import numpy as np;

from Transformer import Constants;

from Transformer.Structure import Structure;


# -----------------
# Sorting Functions
# -----------------

def GroupStructuresBySpacegroup(structures, degeneracies = None, tolerance = None):
    spacegroupGroups = { };

    for i, structure in enumerate(structures):
        # spacegroupGroups will be keyed by tules of (spacegroupNumber, spacegroupSymbol), as returned by the Structure.GetSpacegroup() routine.

        key = structure.GetSpacegroup(tolerance = tolerance);

        if key in spacegroupGroups:
            groupStructures, groupDegeneracies = spacegroupGroups[key];

            groupStructures.append(structure);

            if degeneracies != None:
                groupDegeneracies.append(degeneracies[i]);
        else:
            spacegroupGroups[key] = (
                [structure], [degeneracies[i]] if degeneracies != None else None
                );

    return spacegroupGroups;


# ----------------------------------
# Symmetry "Deduplication" Functions
# ----------------------------------

def MergeStructureSet(structures, degeneracies = None, parentSymmetryOperations = None, tolerance = None, compareLatticeVectors = True, compareAtomTypes = True):
    # If the user does not set a tolerance, set it to the default value used by the Structure class.

    if tolerance == None:
        tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

    # If degeneracies are not supplied, the variable is initialised to a list of ones.
    # In this case, degeneracies will store a count of the number of occurrences of each structure.

    if degeneracies == None:
        degeneracies = [1] * len(structures);

    # Test for and remove duplicates.

    pointer = 0;

    while pointer < len(structures) - 1:
        # Set a reference structure.

        structureRef = structures[pointer];

        # Build a list of structures to compare the remaining ones in the set to.

        compareStructures = None;

        if parentSymmetryOperations != None:
            # If a set of symmetry operations were supplied, apply each one in turn to generate transformed structures to compare to the others in the set.

            compareStructures = [
                structureRef.GetSymmetryTransform(rotation, translation)
                    for rotation, translation in parentSymmetryOperations
                ];
        else:
            # If not, compare to the reference structure itself.

            compareStructures = [structureRef];

        # Loop over structures to compare and remove duplicates.

        for compareStructure in compareStructures:
            # Mark duplicates for removal.

            removeIndices = [];

            for i, structure in enumerate(structures[pointer + 1:]):
                # For testing equality, the atom positions are always compared, and the lattice vectors and atom types may also be compared depending on the parameters.
                # The order of the equivalence tests is set so as to perform the least computationally-demanding ones first.

                remove = False;

                if not remove and compareLatticeVectors:
                    remove = compareStructure.CompareLatticeVectors(structure);

                if not remove and compareAtomTypes:
                    remove = compareStructure.CompareAtomTypeNumbers(structure);

                if not remove:
                    remove = compareStructure.CompareAtomPositions(structure);

                if remove:
                    removeIndices.append(i + pointer + 1);

            # If there are structures to remove, add their degeneracies to that of the reference structure, and remove the entries from the structure set and degeneracy list.

            if len(removeIndices) > 0:
                for index in removeIndices:
                    degeneracies[pointer] = degeneracies[pointer] + degeneracies[index];

                structures = [
                    structure for i, structure in enumerate(structures)
                        if i not in removeIndices
                    ];

                degeneracies = [
                    degeneracy for i, degeneracy in enumerate(degeneracies)
                        if i not in removeIndices
                    ];

        pointer = pointer + 1;

    # Return the merged structure set and associated degeneracies.

    return (structures, degeneracies);


# -----------------
# Utility Functions
# -----------------

def CartesianToFractionalCoordinates(latticeVectors, atomPositions):
    # Treat the lattice vectors as a 3x3 matrix and invert to obtain the transformation matrix to fractional coordinates.

    transformationMatrix = np.linalg.inv(latticeVectors);

    # Return the atom positions multiplied by the transformation matrix.

    return [np.dot(position, transformationMatrix) for position in atomPositions];
