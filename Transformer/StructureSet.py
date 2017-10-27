# Transformer/StructureSet.py


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

class StructureSet(object):
    # -----------
    # Constructor
    # -----------

    def __init__(self, tolerance = None, parentSymmetryOperations = None):
        # If tolerance is not set, use the default from the Structure class.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Store variables.

        self._tolerance = tolerance;
        self._parentSymmetryOperations = parentSymmetryOperations;

        # Initialise an internal dictionary to store structures and degeneracies in spacegroup groups.

        self._structureSet = { };

    # ---------------
    # Private Methods
    # ---------------

    def _AddStructures(self, structures, degeneracies, isUnion = False):
        # Load fields.

        tolerance = self._tolerance;
        parentSymmetryOperations = self._parentSymmetryOperations;

        structureSet = self._structureSet;

        # If isUnion is set, we only need to compare new structures to the ones currently in the set.
        # For this, we need to store the number of structures in each spacegroup group before we add any new ones.

        compareMaxIndices = None;

        if isUnion:
            compareMaxIndices = {
                key : len(structureList) for key, (structureList, _)
                    in structureSet.items()
                };

        # Keep track of the number of new structures added to the set.

        addCount = 0;

        # Loop over new structures and degeneracies.

        for structure, degeneracy in zip(structures, degeneracies):
            # Find the spacegroup.

            spacegroup = structure.GetSpacegroup(tolerance = tolerance);

            if spacegroup not in structureSet:
                # If there is no entry in the structure set for the spacegroup, create one.

                structureSet[spacegroup] = [[structure], [degeneracy]];

                addCount += 1;
            else:
                # Load the lists of structures and degeneracies.

                structureList, degeneracyList = structureSet[spacegroup];

                # Depending on whether compareMaxIndices has been initialised, work out which structures in the list we need to compare against.

                compareMaxIndex = None;

                if compareMaxIndices != None:
                    compareMaxIndex = compareMaxIndices[spacegroup] if spacegroup in compareMaxIndices else 0;
                else:
                    compareMaxIndex = len(structureList);

                # Load the atom count.

                atomCount = structure.GetAtomCount();

                # Varible to keep track of the index of the structure in the set the new structure matches with.

                matchIndex = None;

                # Variable to store symmetry-transformed positions for the new structure.
                # These are (relatively) expensive to generate and in some cases may not be required, so lazy initialisation is used.

                transformedPositions = None;

                # Compare the new structure to those in the current set (up to compareMaxIndex).

                for i in range(0, compareMaxIndex):
                    compareStructure = structureList[i];

                    # Compare atom counts.

                    match = compareStructure.GetAtomCount() == atomCount;

                    # Compare lattice vectors if required.

                    if match:
                        match = structure.CompareLatticeVectors(compareStructure, tolerance = tolerance);

                    # Compare atom-type numbers if required.

                    if match:
                        atomTypeNumbers1 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);
                        atomTypeNumbers2 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);

                        match = (atomTypeNumbers1 == atomTypeNumbers2).all();

                    # Compare atom positions if required.

                    if match:
                        # Initialise transformedPositions if required.

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
                        # For the compareAtomIndexRanges (fourth) parameter of both functions, we set a single range spanning all the atoms.

                        if _Cython:
                            match = _StructureSet._CompareAtomPositions(
                                compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance, [(0, atomCount)]
                                );
                        else:
                            match = _CompareAtomPositions(
                                compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance, [(0, atomCount)]
                                );

                    if match:
                        # If the new structure matches, record the index of the match and break.

                        matchIndex = i;
                        break;

                if matchIndex != None:
                    # If a match was found, update the degeneracy of the matching structure in the set.

                    degeneracyList[matchIndex] += degeneracy;
                else:
                    # If not, add the new structure and degeneracy to the internal sets.

                    structureList.append(structure);
                    degeneracyList.append(degeneracy);

                    addCount += 1;

        # Return the number of structures added to the set.

        return addCount;

    # --------------
    # Public Methods
    # --------------

    def GetStructureSet(self):
        return self._structureSet;

    def GetStructureSetFlat(self):
        # Convert the structure set dictionary into "flat" lists of structures and associated degeneracies.

        structureSet = self._structureSet;

        structuresFlat, degeneraciesFlat = [], [];

        for key in sorted(structureSet.keys()):
            structureList, degeneracyList = structureSet[key]

            structuresFlat = structuresFlat + structureList;
            degeneraciesFlat = degeneraciesFlat + degeneracyList;

        return (structuresFlat, degeneraciesFlat);

    def GetStructureCount(self):
        structureCount = sum(
            len(structureList) for structureList, _ in self._structureSet.values()
            );

        return structureCount;

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

        # Technically, a union operation is only valid when the calling and argument structure sets are configured with the same equality comparisons.
        # If this is the case, we can apply a performance optimisation: we only need to compare structures in the new structure set to those in the calling one _before_ the union, i.e. as we add new structures to the calling set, we do not need to compare them to additional ones added from the new structure set.
        # If not, issue a warning and default to the behaviour Update() method.

        isUnion = self.CompareEquivalenceSettings(structureSet);

        if not isUnion:
            warnings.warn("UpdateUnion() is only valid when the supplied StructureSet object is set to use the same equality-comparison settings as the calling object.", UserWarning);

        return self._AddStructures(*structureSet.GetStructureSetFlat(), isUnion = isUnion);

    def CompareEquivalenceSettings(self, structureSet):
        # Compare tolerances.

        toleranceEquivalenceThreshold = min(self._tolerance, structureSet._tolerance) / 10.0;

        equivalent = (math.fabs(self._tolerance - structureSet._tolerance) < toleranceEquivalenceThreshold);

        if equivalent:
            # Compare parent symmetry operations.

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

        return equivalent;


# ---------------------
# StructureSetOpt Class
# ---------------------

class StructureSetOpt(StructureSet):
    # -----------
    # Constructor
    # -----------

    def __init__(self, expectedAtomCount, tolerance = None, parentSymmetryOperations = None, compareAtomIndexRanges = None):
        # Call the base-class constructor with the tolerance and parent symmetry operations.

        super(StructureSetOpt, self).__init__(tolerance, parentSymmetryOperations);

        # Validate expectedAtomCount.

        if expectedAtomCount <= 0:
            raise Exception("Error: expectedAtomCount must be > 0.");

        # If compareAtomIndexRanges is supplied, make sure the index ranges are sorted and validate.

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

            _, index2 = compareAtomIndexRanges[-1];

            if index2 > expectedAtomCount:
                raise Exception("Error: If supplied, ranges in compareAtomIndexRanges must be consistent with expectedAtomCount.");

        # Store fields.

        self._expectedAtomCount = expectedAtomCount;

        self._compareAtomIndexRanges = compareAtomIndexRanges;

    # ---------------
    # Private Methods
    # ---------------

    def _AddStructures(self, structures, degeneracies, isUnion = False):
        # Override the _AddStructures() method from the StructureSet class.
        # Unfortunately, there was no (obvious) way to do this without some code duplication.

        # Check the structure atom counts.

        expectedAtomCount = self._expectedAtomCount;

        for structure in structures:
            if structure.GetAtomCount() != expectedAtomCount:
                raise Exception("Error: One or more supplied structures do not have the expected number of atoms.");

        # Load fields.

        tolerance = self._tolerance;
        parentSymmetryOperations = self._parentSymmetryOperations;
        compareAtomIndexRanges = self._compareAtomIndexRanges;

        structureSet = self._structureSet;

        # The next chunk of initialisation code is duplicated from the _AddStructure() method from the StructureSet class.

        compareMaxIndices = None;

        if isUnion:
            compareMaxIndices = {
                key : len(structureList) for key, (structureList, _)
                    in structureSet.items()
                };

        # Count the number of added structures.

        addCount = 0;

        # Loop over new structures and degeneracies.

        for structure, degeneracy in zip(structures, degeneracies):
            spacegroup = structure.GetSpacegroup(tolerance = tolerance);

            if spacegroup not in structureSet:
                structureSet[spacegroup] = [[structure], [degeneracy]];

                addCount += 1;
            else:
                structureList, degeneracyList = structureSet[spacegroup];

                # Set compareMaxIndex.

                compareMaxIndex = None;

                if compareMaxIndices != None:
                    compareMaxIndex = compareMaxIndices[spacegroup] if spacegroup in compareMaxIndices else 0;
                else:
                    compareMaxIndex = len(structureList);

                # We will always need the transformed positions, so we don't need to use lazy initialisation.

                transformedPositions = None;

                if self._parentSymmetryOperations != None:
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

                # Index of first matching structure in the set.

                matchIndex = None;

                for i in range(0, compareMaxIndex):
                    compareStructure = structureList[i];

                    # Compare positions.

                    match = None;

                    if _Cython:
                        match = _StructureSet._CompareAtomPositions(
                            compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance, compareAtomIndexRanges
                            );
                    else:
                        match = _CompareAtomPositions(
                            compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance, compareAtomIndexRanges
                            );

                    if match:
                        matchIndex = i;
                        break;

                # If a match was found, update the degeneracy of the matching structure in the set; if not, add the new structure/degeneracy to the lists.

                if matchIndex != None:
                    degeneracyList[matchIndex] += degeneracy;
                else:
                    structureList.append(structure);
                    degeneracyList.append(degeneracy);

                    addCount += 1;

        # Return the number of structures added to the set.

        return addCount;

    # --------------
    # Public Methods
    # --------------

    def CompareEquivalenceSettings(self, structureSetOpt):
        # If the supplied structure set is not of type StructureSetOpt, the number of atoms isn't guarenteed to be fixed -> return False.

        if not isinstance(structureSetOpt, StructureSetOpt):
            return False;

        # Compare atom index ranges, if set.

        if self._compareAtomIndexRanges != None:
            if structureSetOpt._compareAtomIndexRanges != None:
                # The index ranges should be in sort order.

                for (index11, index12), (index21, index22) in zip(self._compareAtomIndexRanges, structureSetOpt._compareAtomIndexRanges):
                    if not (index11 == index21) and (index12 == index22):
                        return False;
            else:
                return False;
        else:
            if structureSetOpt._compareAtomIndexRanges != None:
                return False;

        # Call the base-class CompareEquivalenceSettings() method.

        return super(StructureSetOpt, self).CompareEquivalenceSettings(structureSetOpt);


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
