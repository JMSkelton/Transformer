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

    # Test for and mark duplicates.

    pointer = 0;

    removeIndices = set();

    while pointer < len(structures) - 1:
        if pointer not in removeIndices:
            # Set a reference structure.

            structureRef = structures[pointer];

            # List of symmetry-transformed structures to compare the remaining ones in the set to.
            # We use lazy initialisation to avoid unnecessary work.

            compareStructures = None;

            # Loop over structures to compare and remove duplicates.

            for i in range(pointer + 1, len(structures)):
                if i not in removeIndices:
                    structure = structures[i];

                    if compareStructures == None:
                        # Initialise compareStructures.

                        if parentSymmetryOperations != None:
                            # If a set of symmetry operations were supplied, apply each one in turn to generate transformed structures to compare to the others in the set.

                            compareStructures = [
                                structureRef.GetSymmetryTransform(rotation, translation)
                                    for rotation, translation in parentSymmetryOperations
                                ];

                            # Prune compareStructures to remove duplicates; this can be a big performance boost for high-symmetry parent structures.

                            compareStructures, _ = MergeStructureSet(
                                compareStructures, tolerance = tolerance, compareLatticeVectors = False, compareAtomTypes = False
                                );
                        else:
                            # If not, compare to the reference structure itself.

                            compareStructures = [structureRef];

                    for compareStructure in compareStructures:
                        # For testing equality, the atom positions are always compared, and the lattice vectors and atom types may also be compared depending on the parameters.
                        # The order of the equivalence tests is set so as to perform the least computationally-demanding ones first.

                        remove = True;

                        if compareLatticeVectors:
                            remove = compareStructure.CompareLatticeVectors(structure);

                        if remove and compareAtomTypes:
                            remove = compareStructure.CompareAtomTypeNumbers(structure);

                        if remove:
                            remove = compareStructure.CompareAtomPositions(structure);

                        if remove:
                            # If a structure is marked for removal, its degeneracy is added to that of the reference.

                            degeneracies[pointer] += degeneracies[i];

                            removeIndices.add(i);

                            break;

        pointer += 1;

    # Remove marked structures and degeneracies from the lists.

    if len(removeIndices) > 0:
        structures = [
            structure for i, structure in enumerate(structures)
                if i not in removeIndices
            ];

        degeneracies = [
            degeneracy for i, degeneracy in enumerate(degeneracies)
                if i not in removeIndices
            ];

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
