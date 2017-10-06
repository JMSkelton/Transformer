# Transformer/StructureSet.py by J. M. Skelton


# -------
# Imports
# -------

import math;
import warnings;

import numpy as np;

from Transformer.Structure import Structure;

# Try to import the Cython-optimised _StructureSet module.
# This module provides core inner routines that replace the Python/NumPy implementations in this module.

_Cython = False;

try:
    import pyximport; pyximport.install(setup_args = { 'include_dirs' : np.get_include() });

    from Transformer import _StructureSet;

    _Cython = True;
except ImportError:
    warnings.warn("Optimised merging functions require the pyximport module from the Cython package.", RuntimeWarning);


# ------------------
# StructureSet Class
# ------------------

class StructureSet:
    # -----------
    # Constructor
    # -----------

    def __init__(
            self,
            compareLatticeVectors = True, compareAtomTypeNumbers = True, compareAtomPositions = True,
            tolerance = None, parentSymmetryOperations = None, compareAtomIndexRanges = None
            ):

        # If tolerance is not set, use the default from the Structure class.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # If compareAtomPositionIndexRanges is supplied, make sure the index ranges are sorted and validate.

        if compareAtomIndexRanges != None:
            compareAtomIndexRanges = sorted(
                [(min(index1, index2), max(index1, index2)) for index1, index2 in compareAtomIndexRanges]
                );

            index1Ref, index2Ref = compareAtomIndexRanges[0];

            if index1Ref < 0:
                raise Exception("Error: If supplied, indices in compareAtomIndexRanges must be zero-based (negative indices are not allowed).");

            for index1, index2 in compareAtomIndexRanges[1:]:
                if index1 < index2Ref:
                    raise Exception("Error: If supplied, index pairs in compareAtomIndexRanges must not overlap.");

        # Store variables.

        self._compareLatticeVectors = compareLatticeVectors;
        self._compareAtomTypeNumbers = compareAtomTypeNumbers;
        self._compareAtomPositions = compareAtomPositions;

        self._tolerance = tolerance;
        self._parentSymmetryOperations = parentSymmetryOperations;
        self._compareAtomIndexRanges = compareAtomIndexRanges;

        # Initialise an internal list of structures and degeneracies.

        self._structureSet, self._degeneracySet = [], [];

    # ---------------
    # Private Methods
    # ---------------

    def _AddStructures(self, structures, degeneracies, compareMaxIndex = None):
        # Load variables.

        compareLatticeVectors = self._compareLatticeVectors;
        compareAtomTypeNumbers = self._compareAtomTypeNumbers;
        compareAtomPositions = self._compareAtomPositions;

        tolerance = self._tolerance;
        parentSymmetryOperations = self._parentSymmetryOperations;
        compareAtomIndexRanges = self._compareAtomIndexRanges;

        structureSet, degeneracySet = self._structureSet, self._degeneracySet;

        # If compareAtomIndexRanges is set, get the maximum index for bounds checking.

        maxCompareAtomIndex = None;

        if compareAtomIndexRanges != None:
            maxCompareAtomIndex = max(
                index2 for _, index2 in compareAtomIndexRanges
                );

        # Keep track of the number of new structures added to the set.

        addCount = 0;

        # Loop over new structures and degeneracies.

        for structure, degeneracy in zip(structures, degeneracies):
            # Load the atom count.

            atomCount = structure.GetAtomCount();

            # If compareAtomPositionIndexRanges has not been set, set it to compare all atom positions, based on the number of atoms in the new structure.
            # If it has, check the maximum index in the ranges does not exceed the number of atoms in the structure to be added.

            if compareAtomIndexRanges == None:
                compareAtomIndexRanges = [(0, atomCount)];
            else:
                if maxCompareAtomIndex != None and maxCompareAtomIndex > atomCount:
                    raise Exception("Error: One or more of the index ranges in compareAtomIndexRanges is out of bounds for an added structure.");

            # Varible to keep track of the index of the structure in the set the new structure matches with.

            matchIndex = None;

            # Variable to store symmetry-transformed positions for the new structure.
            # These are (relatively) expensive to generate and in some cases may not be required, so lazy initialisation is used.

            transformedPositions = None;

            # Compare the new structure to those in the current set.
            # If compareMaxIndex is set, only compare the new structure against the first compareMaxIndex structures in the current set.

            for i, compareStructure in enumerate(structureSet if compareMaxIndex == None else structureSet[:compareMaxIndex]):
                # If two structures do not contain the same number of atoms, they are automatically identified as being different.

                match = compareStructure.GetAtomCount() == structure.GetAtomCount();

                # Compare lattice vectors if required.

                if match and compareLatticeVectors:
                    match = structure.CompareLatticeVectors(compareStructure, tolerance = tolerance);

                # Compare atom type numbers if required.

                if match and compareAtomTypeNumbers:
                    atomTypeNumbers1 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);
                    atomTypeNumbers2 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);

                    for index1, index2 in compareAtomIndexRanges:
                        match = (atomTypeNumbers1[index1:index2] == atomTypeNumbers2[index1:index2]).all();

                        if not match:
                            break;

                # Compare atom positions, if required.

                if match and compareAtomPositions:
                    # Initialise transformedPositions, if required.

                    if transformedPositions is None:
                        # If a set of symmetry operations for a parent structure have been supplied, generate symmetry-transformed positions for the new structure; if not, compare the positions directly.

                        if self._parentSymmetryOperations != None:
                            # Use the Cython-optimised routines if available; if not, fall back to the NumPy implementation.

                            if _Cython:
                                transformedPositions = _StructureSet._GenerateSymmetryTransformedPositions(
                                    structure, parentSymmetryOperations, tolerance
                                    );
                            else:
                                transformedPositions = _GenerateSymmetryTransformedPositions(
                                    structure, parentSymmetryOperations, tolerance
                                    );
                        else:
                            transformedPositions = np.array(
                                [structure.GetAtomPositionsNumPy(copy = False)]
                                );

                    # Use the Cython-optimised comparison routines if available; if not, fall back to the NumPy implementation.

                    if _Cython:
                        match = _StructureSet._CompareAtomPositions(compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance, compareAtomIndexRanges);
                    else:
                        match = _CompareAtomPositions(compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance, compareAtomIndexRanges);

                if match:
                    # If the new structure matches, record the index of the match and break.

                    matchIndex = i;
                    break;

            if matchIndex != None:
                # If a match was found, update the degeneracy of the matching structure in the set.

                degeneracySet[matchIndex] += degeneracy;
            else:
                # If not, add the new structure and degeneracy to the internal sets.

                structureSet.append(structure);
                degeneracySet.append(degeneracy);

                addCount += 1;

        # Return the number of structures added to the set.

        return addCount;

    # --------------
    # Public Methods
    # --------------

    def GetStructures(self):
        return self._structureSet;

    def GetDegeneracies(self):
        return self._degeneracySet;

    def GetStructureCount(self):
        return len(self._structureSet);

    def Add(self, structure, degeneracy = 1):
        # Sanity check.

        if structure == None:
            raise Exception("Error: structure cannot be None.");

        addCount = self._AddStructures([structure], [degeneracy]);

        return True if addCount == 1 else False;

    def Update(self, structures, degeneracies = None):
        # Sanity checks.

        if structures == None:
            raise Exception("Error: structures cannot be None.");

        if degeneracies != None and len(degeneracies) != len(structures):
            raise Exception("Error: If supplied, degeneracies must have the same length as structures.");

        if degeneracies == None:
            degeneracies = [1] * len(structures);

        return self._AddStructures(structures, degeneracies);

    def UpdateUnion(self, structureSet):
        if structureSet == None:
            raise Excpetion("Error: structureSet cannot be None.");

        compareMaxIndex = None;

        # Technically, a union operation is only valid when the calling and argument structure sets are configured with the same equality comparisons.
        # If this is the case, we can apply a performance optimisation: we only need to compare structures in the new structure set to those in the calling one _before_ the union, i.e. as we add new structures to the calling set, we do not need to compare them to additional ones added from the new structure set.
        # If not, issue a warning and default to the behaviour Update() method.

        if self.CompareEquivalenceSettings(structureSet):
            compareMaxIndex = len(self._structureSet);
        else:
            warnings.warn("UpdateUnion() is only valid when the supplied StructureSet object is set to use the same equality-comparison settings as the calling object.", UserWarning);

        return self._AddStructures(structureSet.GetStructures(), structureSet.GetDegeneracies(), compareMaxIndex = compareMaxIndex);

    def CompareEquivalenceSettings(self, structureSet):
        equivalent = True;

        equivalent = equivalent and (self._compareLatticeVectors == structureSet._compareLatticeVectors);
        equivalent = equivalent and (self._compareAtomTypeNumbers == structureSet._compareAtomTypeNumbers);
        equivalent = equivalent and (self._compareAtomPositions == structureSet._compareAtomPositions);

        if equivalent:
            # Compare tolerances.

            toleranceEquivalenceThreshold = min(self._tolerance, structureSet._tolerance) / 10.0;

            equivalent = equivalent and (math.fabs(self._tolerance - structureSet._tolerance) < toleranceEquivalenceThreshold);

        if equivalent:
            # Compare parent symmetry operations, if set.

            if self._parentSymmetryOperations != None:
                if structureSet._parentSymmetryOperations != None:
                    # Use the set tolerance (should be the same for the calling and argument StructureSet objects) to compare translations.

                    tolerance = self._tolerance;

                    # Compare symmetry operations.

                    for (rotation1, translation1), (rotation2, translation2) in zip(self._parentSymmetryOperations, structureSet._parentSymmetryOperations):
                        equivalent = equivalent and (rotation1 == rotation2).all();
                        equivalent = equivalent and (np.abs(translation1 - translation2) < tolerance).all();

                        if not equivalent:
                            break;
                else:
                    equivalent = False;
            else:
                if structureSet._parentSymmetryOperations != None:
                    equivalent = False;

        if equivalent:
            # Compare atom index ranges, if set.

            if self._compareAtomIndexRanges != None:
                if structureSet._compareAtomIndexRanges != None:
                    # The index ranges should be in sort order.

                    for (index11, index12), (index21, index22) in zip(self._compareAtomIndexRanges, structureSet._compareAtomIndexRanges):
                        equivalent = equivalent and (index11 == index21);
                        equivalent = equivalent and (index12 == index22);

                        if not equivalent:
                            break;
                else:
                    equivalent = False;
            else:
                if structureSet._compareAtomIndexRanges != None:
                    equivalent = False;

        return equivalent;


# --------------
# Static Methods
# --------------

# These should ideally be included in the StructureSet class, but Python 2.x doesn't allow classes to have static methods attached to them.

def _CompareAtomPositions(comparePositions, refTransformedPositions, tolerance, compareAtomIndexRanges):
    # Start by assuming all the sets of transformed positions match.

    match = np.ones(refTransformedPositions.shape[0], dtype = np.bool);

    # Loop over index ranges in compareAtomIndexRanges.

    for index1, index2 in compareAtomIndexRanges:
        # For each, perform a logical and with match.

        compareResult = np.all(
            np.abs(refTransformedPositions[:, index1:index2, :] - comparePositions[index1:index2, :]) < tolerance, axis = (1, 2)
            );

        match = np.logical_and(match, compareResult);

    # If any elements in match remain True, comparePositions matches one or more elements in refTransformedPositions.

    return match.any();

def _GenerateSymmetryTransformedPositions(structure, symmetryOperations, tolerance):
    # Get the atom data from the supplied Structure object.

    atomData = structure.GetAtomDataNumPy(copy = False);

    numAtoms = len(atomData);
    numSymOps = len(symmetryOperations);

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

    # In rare cases, applying symmetry operations yields fractional coordinates very close to 1.0 that aren't wrapped to 0.0 by the % operator.
    # When comparing structures, this can cause equivalent sets of positions to be erroneously identified as different.
    # To work around this, we explicitly set coordinates within the symmetry tolerance of 1.0 to 0.0.

    transformedPositions[np.abs(transformedPositions - 1.0) < tolerance] = 0.0;

    # Sort the data blocks.

    for i in range(0, numSymOps):
        transformedStructures[i].sort();

    # Eliminate redundant transformed structures.
    # This is a relatively expensive operation, but appears to be worth it when using the pure-Python merging routines.

    transformedStructures = np.unique(transformedStructures, axis = 0)

    # Return a view to the transformed positions.

    return transformedStructures.view(dtype = np.float64).reshape((len(transformedStructures), structure.GetAtomCount(), 4))[:, :, 1:];
