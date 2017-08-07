# Transformer/StructureTools.py by J. M. Skelton


# -------
# Imports
# -------

import warnings;

import numpy as np;

from Transformer import Constants;

from Transformer.Structure import Structure;

# Try to import the tqdm function to display a progress bar in the MergeStructureSet() routine.

ImportedTQDM = False;

try:
    from tqdm import tqdm;

    ImportedTQDM = True;
except ImportError:
    warnings.warn("The tqdm module could not be imported -> displaying progress bars in MergeStructureSet() will be disabled.", RuntimeWarning);

    pass;


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

def MergeStructureSet(structures, degeneracies = None, parentSymmetryOperations = None, tolerance = None, compareLatticeVectors = True, compareAtomTypes = True, progressBar = False):
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

    # Set up a primary iterator.
    # If progressBar is set and the tqdm module is available, wrap it in the tqdm() function to display a progress bar.

    iValues = range(0, len(structures));

    if ImportedTQDM and progressBar:
        iValues = tqdm(iValues);

    # Loop over reference structures.

    for i in iValues:
        if i not in removeIndices:
            structureRef = structures[i];

            # List of symmetry-transformed structures to compare the remaining ones in the set to.
            # We use lazy initialisation to avoid unnecessary work.

            comparePositions = None;

            # Loop over structures to compare and remove duplicates.

            # Testing shows that applying the symmetry operations of a parent structure (at least as implemented here) can be "assymetric", i.e. if trans(A) == B is false, trans(B) == A might be true.
            # Therefore, for each reference structure, we compare to _all_ other structures that haven't been eliminated yet.

            for j in range(0, len(structures)):
                if j != i and j not in removeIndices:
                    structure = structures[j];

                    if comparePositions is None:
                        # Initialise compareStructures.

                        transformedStructures = None;

                        if parentSymmetryOperations != None:
                            # If a set of symmetry operations have been supplied, generate a set of symmetry-transformed structures to compare to.

                            transformedStructures = GenerateSymmetryTransformedStructures(
                                structureRef, parentSymmetryOperations
                                );
                        else:
                            # If not, compare against the original reference structure.

                            transformedStructures = np.array(
                                [structureRef.GetAtomDataNumPy()]
                                );

                        comparePositions = transformedStructures.view(dtype = np.float64).reshape((len(transformedStructures), structureRef.GetAtomCount(), 4))[:, :, 1:];

                    # If the compareLatticeVectors and/or compareAtomTypes flags are set, compare the lattice vectors/atom-type numbers first.
                    # Even if the reference structure was "expanded" by symmetry operations, these comparisons onlt need to be performed once.

                    remove = True;

                    if compareLatticeVectors:
                        remove = structureRef.CompareLatticeVectors(structure);

                    if remove and compareAtomTypes:
                        remove = structureRef.CompareAtomTypeNumbers(structure);

                    if remove:
                        compareResult = np.all(
                            np.abs(comparePositions - structure.GetAtomPositionsNumPy(copy = False)) < tolerance, axis = (1, 2)
                            );

                        remove = compareResult.any();

                    if remove:
                        # If a structure is marked for removal, its degeneracy is added to that of the reference.

                        degeneracies[i] += degeneracies[j];

                        removeIndices.add(j);

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

def GenerateSymmetryTransformedStructures(structure, symmetryOperations):
    # Get the atom data from the supplied Structure object.

    atomData = structure.GetAtomDataNumPy(copy = False);

    numAtoms = len(atomData);
    numSymOps = len(symmetryOperations);

    # Create an N_ops x N_a NumPy array of the Structure._AtomDataType data structure and initialise it with the atom data.

    transformedStructures = np.zeros(
        (numSymOps, numAtoms), dtype = Structure._AtomDataType
        );

    transformedStructures[:] = atomData;

    transformedPositions = transformedStructures.view(dtype = np.float64).reshape((numSymOps, numAtoms, 4))[:, :, 1:];

    # Loop over symmetry operations and transform the positions in each block of data.

    for i, (rotation, translation) in enumerate(symmetryOperations):
        # Apply the rotation and translation.

        transformedPositions[i] = np.einsum('jk,ik', rotation, transformedPositions[i]) + translation;

    # Clamp the modified positions to the range [0, 1].

    transformedPositions %= 1.0;

    # Sort the data blocks.

    for i in range(0, numSymOps):
        transformedStructures[i].sort();

    # Eliminate redundant transformed structures.
    # With small numbers of symmetry operations, this provides little benefit, but does not seem to degrade performance either.

    return np.unique(transformedStructures, axis = 0);
