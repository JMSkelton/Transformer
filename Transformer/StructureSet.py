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
    warnings.warn("Optimised merging requires the pyximport module from the Cython package.", RuntimeWarning);


# ------------------
# StructureSet Class
# ------------------

class StructureSet(object):
    """
    Class for maintaining a set of unique Structure objects and associated degeneracies.

    Structures can be compared with a set of user-supplied symmetry operations, which are used to transform structures to symmetry-equivalent configurations.

    Internally, the set is stored in a dictionary where:
        - the keys are (spacegroup_number, spacegroup_symbol) tuples; and
        - the values are tuples of (structures, degeneracies) lists

    If the pyximport module from the Cython package is available, the symmetry transformations and structure comparisons are sped up using C kernels.
    """

    # -----------
    # Constructor
    # -----------

    def __init__(
        self,
        tolerance = None, parentSymmetryOperations = None,
        structures = None, degeneracies = None, noInitialMerge = False
        ):

        """
        Class constructor.

        Keyword arguments:
            tolerance -- symmetry tolerance for structure comparisons.
            parentSymmetryOperations -- symmetry operations to be used to transform structures to symmetry-equivalent configurations.
            structures, degeneracies -- lists of structures and degeneracies to initialise the set with.
            noInitialMerge -- if True, do not merge the initial structures (if supplied).
        """

        # If tolerance is not set, use the default from the Structure class.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # If an iniital list of structures is provided, check/initialise the list of degeneracies.

        if structures != None:
            if degeneracies != None:
                if len(degeneracies) != len(structures):
                    raise Exception("Error: If supplied, degeneracies must contain one entry for each structure in the list of structures.");
            else:
                # If degeneracies are not supplied, set it to a list of ones.

                degeneracies = [1] * len(structures);
        else:
            if degeneracies != None:
                raise Exception("Error: degeneracies can only be supplied alongside a list of structures.");

        self._tolerance = tolerance;
        self._parentSymmetryOperations = parentSymmetryOperations;

        self._structureSet = { };

        # If a list of structure has been supplied, initialise the structure set.

        if structures != None:
            if noInitialMerge:
                # If noInitialMerge is set, bypass merging and populate the set directly.

                structureSet = self._structureSet;

                for structure, degeneracy in zip(structures, degeneracies):
                    spacegroup = structure.GetSpacegroup();

                    if spacegroup not in structureSet:
                        structureSet[spacegroup] = ([structure], [degeneracy]);
                    else:
                        structureList, degeneracyList = structureSet[spacegroup];

                        structureList.append(structure);
                        degeneracyList.append(degeneracy);
            else:
                # If not, pass the structures and degeneracies to the _AddStructures() method to be merged in.

                self._AddStructures(structures, degeneracies, isUnion = False);


    # ---------------
    # Private Methods
    # ---------------

    def _AddStructures(self, structures, degeneracies, isUnion):
        """
        Merge a list of structures and degeneracies into the set and return the number of structures added.
        If isUnion is True, assume the structures in the list are unique, and only compare them to structures in the initial set while merging.
        """

        tolerance = self._tolerance;
        parentSymmetryOperations = self._parentSymmetryOperations;

        structureSet = self._structureSet;

        # If isUnion is set, we assume the new structures are unique and apply a performance optimisation by only comparing new structures to those in the initial set.
        # For this, we need to store the number of structures in each spacegroup group before we add any new ones.

        compareMaxIndices = None;

        if isUnion:
            compareMaxIndices = {
                key : len(structureList) for key, (structureList, _)
                    in structureSet.items()
                };

        # Keep track of the number of new structures added to the set.

        addCount = 0;

        for structure, degeneracy in zip(structures, degeneracies):
            atomCount = structure.GetAtomCount();
            spacegroup = structure.GetSpacegroup(tolerance = tolerance);

            if spacegroup not in structureSet:
                # If there is no entry in the structure set for the spacegroup, create one.

                structureSet[spacegroup] = ([structure], [degeneracy]);

                addCount += 1;

            else:
                structureList, degeneracyList = structureSet[spacegroup];

                # Work out the range of structures to compare against.

                compareMaxIndex = None;

                if compareMaxIndices != None:
                    compareMaxIndex = compareMaxIndices[spacegroup] if spacegroup in compareMaxIndices else 0;
                else:
                    compareMaxIndex = len(structureList);

                # Store the index of the existing structure the new structure matches.

                matchIndex = None;

                # Variable to store symmetry-transformed positions for the new structure.
                # These are (relatively) expensive to generate and in some cases may not be required, so lazy initialisation is used.

                transformedPositions = None;

                # Compare against structures in the current set up to compareMaxIndex.

                for i in range(0, compareMaxIndex):
                    compareStructure = structureList[i];

                    # Compare atom counts.

                    match = compareStructure.GetAtomCount() == atomCount;

                    # Compare lattice vectors.

                    if match:
                        match = structure.CompareLatticeVectors(compareStructure, tolerance = tolerance);

                    # Compare atom-type numbers.

                    if match:
                        atomTypeNumbers1 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);
                        atomTypeNumbers2 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);

                        match = (atomTypeNumbers1 == atomTypeNumbers2).all();

                    # Compare atom positions.

                    if match:
                        if transformedPositions is None:
                            if self._parentSymmetryOperations != None:
                                # If a set of symmetry operations have been supplied, generate symmetry-transformed positions for the new structure.

                                transformedPositions = _StructureSet._GenerateSymmetryTransformedPositions(
                                    structure, parentSymmetryOperations, tolerance
                                    );
                            else:
                                # If not, just compare positions directly.

                                transformedPositions = np.array(
                                    [structure.GetAtomPositionsNumPy(copy = False)]
                                    );

                        # For the compareAtomIndexRanges (fourth) parameter, we set a single range spanning all the atoms.

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
                    # If not, add the new structure and degeneracy to the lists.

                    structureList.append(structure);
                    degeneracyList.append(degeneracy);

                    addCount += 1;

        # Return the number of structures added to the set.

        return addCount;

    # --------------
    # Public Methods
    # --------------

    def GetStructureSet(self):
        """
        Return the internal dictionary representation of the structure set.
        The keys are (spacegroup_number, spacegroup_symbol) tuples, and the values are tuples of (structures, degeneracies) lists.
        """

        return self._structureSet;

    def GetStructureSetFlat(self):
        """ Return the internal structure set as "flat" lists of structures and degeneracies. """

        structureSet = self._structureSet;

        structuresFlat, degeneraciesFlat = [], [];

        for key in sorted(structureSet.keys()):
            structureList, degeneracyList = structureSet[key]

            structuresFlat = structuresFlat + structureList;
            degeneraciesFlat = degeneraciesFlat + degeneracyList;

        return (structuresFlat, degeneraciesFlat);

    def GetStructureCount(self):
        """ Return the number of structures in the set. """

        structureCount = sum(
            len(structureList) for structureList, _ in self._structureSet.values()
            );

        return structureCount;

    def Add(self, structure, degeneracy = 1):
        """
        Merge structure into the set with the optional supplied degeneracy.
        Returns True if the structure was added as a new structure, and False if it was merged with one already in the set.
        """

        if structure == None:
            raise Exception("Error: structure cannot be None.");

        addCount = self._AddStructures([structure], [degeneracy], isUnion = False);

        return addCount == 1;

    def Update(self, structures, degeneracies = None):
        """ Merge the list of structures into the set with the optional supplied degeneracies, and return the number of new structures added. """

        if structures == None:
            raise Exception("Error: structures cannot be None.");

        if degeneracies != None and len(degeneracies) != len(structures):
            raise Exception("Error: If supplied, degeneracies must have the same length as structures.");

        if degeneracies == None:
            degeneracies = [1] * len(structures);

        return self._AddStructures(structures, degeneracies, isUnion = False);

    def UpdateUnion(self, structureSet):
        """
        Perform a union operation with structureSet and return the number of new structures added.
        A union is only valid when both structure sets use the same parameters for identifying unique structures.
        """

        if structureSet == None:
            raise Excpetion("Error: structureSet cannot be None.");

        # Check the union is valid and issue a warning if not.

        isUnion = self.CompareEquivalenceSettings(structureSet);

        if not isUnion:
            warnings.warn("UpdateUnion() is only valid when the supplied structure set is set to use the same comparison parameters as the calling one.", UserWarning);

        return self._AddStructures(*structureSet.GetStructureSetFlat(), isUnion = isUnion);

    def CompareEquivalenceSettings(self, structureSet):
        """ Compare the parameters used to determine structural equality with those of structureSet. """

        # Compare tolerances.

        tolerance = self._tolerance;

        # To compare tolerances, we set a threshold an order of magnitude tighter than the minimum tolerance.

        equivalent = math.fabs(tolerance - structureSet._tolerance) < (min(tolerance, structureSet._tolerance) / 10.0);

        if equivalent:
            # Compare symmetry operations.

            if self._parentSymmetryOperations != None:
                if structureSet._parentSymmetryOperations != None:
                    for (rotation1, translation1), (rotation2, translation2) in zip(self._parentSymmetryOperations, structureSet._parentSymmetryOperations):
                        # rotation matrices are integers, and translation vectors are floating-point numbers.

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

    def CloneNew(self, structures = None, degeneracies = None, noInitialMerge = False):
        """
        Return a new StructureSet object with the same comparison settings.
        Optionally, pass the structures, degeneracies and noInitialMerge argument to the constructor.
        """

        return StructureSet(
            tolerance = self._tolerance, parentSymmetryOperations = self._parentSymmetryOperations,
            structures = structures, degeneracies = degeneracies, noInitialMerge = noInitialMerge
            );


# ---------------------
# StructureSetOpt Class
# ---------------------

class StructureSetOpt(StructureSet):
    """
    An optimised StructureSet class which:
        - assumes all the structures have the same lattice vectors and composition; and
        - can optionally be configured to compare only certain ranges of atom positions

    These performance optimisations are used by the atomic-substitutions framework routines.
    """

    # -----------
    # Constructor
    # -----------

    def __init__(
        self,
        expectedAtomCount, tolerance = None, parentSymmetryOperations = None,
        structures = None, degeneracies = None, noInitialMerge = False,
        compareAtomIndexRanges = None
        ):

        """
        Constructor.

        Arguments:
            expectedAtomCount -- number of atoms in structures.

        Keyword arguments:
            tolerance, parentSymmetryOperations -- as for StructureSet class.
            structures, degeneracies, noInitialMerge -- as for StructureSet class.
            compareAtomIndexRanges -- tuples of (start_inclusive, end_exclusive) index ranges over which to compare atoms.
        """

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

        self._expectedAtomCount = expectedAtomCount;
        self._compareAtomIndexRanges = compareAtomIndexRanges;

        # Call the base-class constructor.

        super(StructureSetOpt, self).__init__(tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, structures = structures, degeneracies = degeneracies, noInitialMerge = noInitialMerge);

    # ---------------
    # Private Methods
    # ---------------

    def _AddStructures(self, structures, degeneracies, isUnion):
        """ Override of the _AddStructures() method from the StructureSet class. """

        # Check atom counts.

        expectedAtomCount = self._expectedAtomCount;

        for structure in structures:
            if structure.GetAtomCount() != expectedAtomCount:
                raise Exception("Error: One or more supplied structures do not have the expected number of atoms.");

        tolerance = self._tolerance;
        parentSymmetryOperations = self._parentSymmetryOperations;
        compareAtomIndexRanges = self._compareAtomIndexRanges;

        structureSet = self._structureSet;

        # Initialisation code duplicated from the _AddStructure() method from the StructureSet class.

        compareMaxIndices = None;

        if isUnion:
            compareMaxIndices = {
                key : len(structureList) for key, (structureList, _)
                    in structureSet.items()
                };

        # Keep track of the number of structures added.

        addCount = 0;

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

                # We will always need the transformed positions -> no need to use lazy initialisation.

                transformedPositions = None;

                if self._parentSymmetryOperations != None:
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
        """
        Compare the parameters used to determine structural equality with those of structureSetOpt.

        Notes:
            This overrides the CompareEquivalenceSettings() method of the StructureSet class, and will return False if structureSetOpt is not an instance of StructureSetOpt.
        """

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

    def CloneNew(self, structures = None, degeneracies = None, noInitialMerge = False, compareAtomIndexRanges = None):
        """
        Return a new StructureSetOpt object with the same comparison settings.
        Optionally, pass the structures, degeneracies, noInitialMerge and compareAtomIndexRanges arguments to the constructor.

        Notes:
            This overrides the CloneNew() method of the StructureSet class.
            The atom-index ranges are not cloned by default, since these are intended as a performance optimisation rather than a comparison setting.
        """

        return StructureSetOpt(
            self._expectedAtomCount, tolerance = self._tolerance, parentSymmetryOperations = self._parentSymmetryOperations,
            structures = structures, degeneracies = degeneracies, noInitialMerge = noInitialMerge,
            compareAtomIndexRanges = compareAtomIndexRanges
            );

# ------------------
# Internal Functions
# ------------------

def _CompareAtomPositions(comparePositions, refTransformedPositions, tolerance, compareAtomIndexRanges):
    """
    Return true if the atomic positions in comparePositions match any of the set of reference symmetry-transformed positions in refTransformedPositions to within the supplied tolerane.
    Only the position ranges is .

    Arguments:
        comparePositions -- an Nx3 NumPy matrix containing a set of atom positions.
        refTransformedPositions -- an MxNx3 NumPy matrix containing a set of reference symmetry-transformed atom positions.
        tolerance -- symmetry tolerance for comparing atom positions.
        compareAtomIndexRanges --

    Notes:
        No default initialisation of compareAtomIndexRanges is performed by this function.
        No parameter or bounds checking is performed by this function.
    """

    # If the Cython-optimised routines in the _StructureSet module are available, use those.

    if _Cython:
        return _StructureSet._CompareAtomPositions(comparePositions, refTransformedPositions, tolerance, compareAtomIndexRanges);

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
    """
    Generate M sets of symmetry-transformed positions from structure using the supplied symmetry operations.

    Arguments:
        structure -- a Structure object with the positions to be transformed.
        symmetryOperations -- a list of (rotation, translation) tuples specifying the symmetry transformations to apply.
        tolerance -- symmetry tolerance for adjusting transformed positions to minimise numerical noise.

    Return value:
        An MxNx3 NumPy matrix containing M transformed sets of N atomic positions specified as NumPy vectors.
        The positions are returned sorted by atom type then position, i.e. they match the internal data layout used in the Structure class.
    """

    # If the Cython-optimised routines in the _StructureSet module are available, use those.

    if _Cython:
        return _StructureSet._GenerateSymmetryTransformedPositions(structure, symmetryOperations, tolerance);

    atomData = structure.GetAtomDataNumPy(copy = False);

    numAtoms = len(atomData);
    numSymOps = len(symmetryOperations);

    transformedStructures = np.zeros(
        (numSymOps, numAtoms), dtype = Structure._AtomDataType
        );

    transformedStructures[:] = atomData;

    transformedPositions = transformedStructures.view(dtype = np.float64).reshape((numSymOps, numAtoms, 4))[:, :, 1:];

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

    # Return a view to the transformed positions.

    return transformedStructures.view(dtype = np.float64).reshape((len(transformedStructures), structure.GetAtomCount(), 4))[:, :, 1:];
