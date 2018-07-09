# Transformer/StructureSet.py


# -------
# Imports
# -------

import math;
import warnings;

import numpy as np;

from Transformer.Structure import Structure;
from Transformer.Utilities import MultiprocessingHelper;

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

    Structure equivalence testing is controlled by a set of parameters which allow for:
        - transforming structures to symmetry-equivalent configurations using a supplied set of symmetry operations;
        - setting a tolerance for comparing lattice vectors and/or positions;
        - comparing any or all of the lattice vectors, atom-type numbers and atom positions; and
        - restricting the comparison to certain ranges of atom indices (if all structures can be expected to have the same numbers of atoms).

    (Several of the above parameters are designed as performance optimisations to avoid unnecessary work in specific cases.)

    Internally, the set is stored in a dictionary where:
        - the keys are (spacegroup_number, spacegroup_symbol) tuples; and
        - the values are tuples of (structures, degeneracies) lists.

    If the pyximport module from the Cython package is available, the symmetry transformations and structure comparisons are sped up using C kernels.
    """

    # -----------
    # Constructor
    # -----------

    def __init__(
        self,
        structures = None, degeneracies = None, noInitialMerge = False,
        tolerance = None, symmetryExpansion = 'none', parentSymmetryOperations = None,
        compareLatticeVectors = True, compareAtomTypeNumbers = True, compareAtomPositions = True,
        expectedAtomCount = None, compareAtomIndexRanges = None
        ):

        """
        Class constructor.

        Keyword arguments:
            structures, degeneracies -- lists of structures and degeneracies to initialise the set with.
            noInitialMerge -- if True, do not merge the initial structures (if supplied).

            tolerance -- symmetry tolerance for structure comparisons.
            symmetryExpansion -- control the coverage of the symmetry expansion when comparing structures ('none', 'fast', 'full'; default: 'none').
            parentSymmetryOperations -- in conjunction with symmetryExpansion, sets the symmetry operations to be used to transform structures to symmetry-equivalent configurations.

            compareLatticeVectors -- if True (default), compare lattice vectors when testing structure equivalence.
            compareAtomTypeNumbers -- if True (default), compare atom-type numbers when testing structure equivalence.
            compareAtomPositions -- if True (default), compare atom positions when testing structure equivalence.

            expectedAtomCount -- if supplied, requires that all structures contain the same number of atoms (used in conjunction with compareAtomIndexRanges).
            compareAtomIndexRanges -- if comparing atom-type numbers and/or atom positions, optionally provides a list of (start_inclusive, end_exclusive) tuples specifying ranges of atom indices to compare (requires expectedAtomCount to be set).

        Notes:
            The symmetryExpansion keyword controls how symmetry operations are used to compare structures when adding new ones to the set.

            Possible options are:
                'none' - do not use symmetry;
                'fast' - expand new structures (but not those in the set) using symmetry operations; and
                'full' - expand new structures and structures in the set using symmetry operations.

            If 'fast' or 'full' are selected without supplying a set of symmetry operations, an error is thrown.
            For efficiency, 'full' caches the expansions of structures in the set, which increases memory usage.
            Preliminary testing suggests that the 'fast' setting strikes a good balance between speed and good symmetry reduction for many problems.

            If noInitialMerge is specified when expectedAtomCount is set, the numbers of atoms in the initial set of structures (if supplied) are checked.
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

        symmetryExpansion = symmetryExpansion.lower();

        if symmetryExpansion not in ['none', 'fast', 'medium', 'full']:
            raise Exception("Error: Unrecognised symmetry-expansion option '{0}'.".format(symmetryExpansion));

        if symmetryExpansion != 'none' and parentSymmetryOperations == None:
            raise Exception("Error: If symmetry expansion is enabled (symmetryExpansion != 'none'), a set of parent symmetry operations must be supplied.");

        self._symmetryExpansion = symmetryExpansion;
        self._parentSymmetryOperations = parentSymmetryOperations;

        self._compareLatticeVectors = compareLatticeVectors;
        self._compareAtomTypeNumbers = compareAtomTypeNumbers;
        self._compareAtomPositions = compareAtomPositions;

        if expectedAtomCount != None:
            if expectedAtomCount <= 0:
                raise Exception("Error: If supplied, expectedAtomCount must be > 0.");

                if structures != None:
                    for structure in structures:
                        if structure.GetAtomCount() != expectedAtomCount:
                            raise Exception("Error: One or more supplied initial structures does not have the number of atoms set by expectedAtomCount.");

        if compareAtomIndexRanges != None:
            if expectedAtomCount == None:
                raise Exception("Error: compareAtomIndexRanges can only be specified when expectedAtomCount is set.");

            # Sorting compareAtomIndexRanges before storing makes comparing equivalence settings in the CompareEquivalenceSettings() method more reliable.

            compareAtomIndexRanges = sorted(compareAtomIndexRanges);

            indices = [index1 for index1, _ in compareAtomIndexRanges] + [index2 for _, index2 in compareAtomIndexRanges];

            if min(indices) < 0 or max(indices) > expectedAtomCount:
                raise Exception("Error: Ranges specified in compareAtomIndexRanges must be between 0 and expectedAtomCount.");

        self._expectedAtomCount = expectedAtomCount;
        self._compareAtomIndexRanges = compareAtomIndexRanges;

        self._structureSet = { };

        symmetryExpansionsCache = None;

        if symmetryExpansion == 'full':
            symmetryExpansionsCache = { };

        self._symmetryExpansionsCache = symmetryExpansionsCache;

        # If a list of structure has been supplied, initialise the structure set.

        if structures != None:
            # if noInitialMerge is set, this is equivalent to performing a union with the internal structure set empty.

            self._AddStructures(structures, degeneracies, isUnion = noInitialMerge);


    # ---------------
    # Private Methods
    # ---------------

    def _FindStructure(self, structure):
        """
        Attempt to find a supplied structure in the internal structure set.

        Return value:
            If the structure is found, a (spacegroup, structure_number) tuple indexing the matched structure in the internal structure set; if not, returns None.
        """

        tolerance = self._tolerance;

        symmetryExpansion = self._symmetryExpansion;
        parentSymmetryOperations = self._parentSymmetryOperations;

        compareLatticeVectors = self._compareLatticeVectors;
        compareAtomTypeNumbers = self._compareAtomTypeNumbers;
        compareAtomPositions = self._compareAtomPositions;

        expectedAtomCount = self._expectedAtomCount;
        compareAtomIndexRanges = self._compareAtomIndexRanges;

        structureSet = self._structureSet;

        symmetryExpansionsCache = self._symmetryExpansionsCache;

        # If an expected atom count has been set, check structure has the correct number of atoms.

        atomCount = structure.GetAtomCount();

        if expectedAtomCount != None and atomCount != expectedAtomCount:
            raise Exception("Error: structure does not have the set expected atom count ({0} != {1}).".format(atomCount, expectedAtomCount));

        # Get the spacegroup of structure and check whether the set contains structures with the same spacegroup.

        spacegroup = structure.GetSpacegroup(tolerance = tolerance);

        if spacegroup not in structureSet:
            return None;

        # Compare structure to those in the set with the same spacegroup.

        structureList, _ = structureSet[spacegroup];

        symmetryExpansionsList = None;

        if symmetryExpansionsCache != None:
            symmetryExpansionsList = symmetryExpansionsCache[spacegroup];

        # List index of matching structure.

        matchIndex = None;

        # Symmetry-transformed positions for the new structure.
        # Since these are (relatively) expensive to generate and are not necessarily required, we use lazy initialisation.

        transformedPositions = None;

        # Loop over structures.

        for i, compareStructure in enumerate(structureList):
            match = True;

            # Compare lattice vectors if required.

            if compareLatticeVectors:
                latticeVectors1 = structure.GetLatticeVectorsNumPy(copy = False);
                latticeVectors2 = compareStructure.GetLatticeVectorsNumPy(copy = False);

                match = np.all(np.abs(latticeVectors1 - latticeVectors2) < tolerance);

            # Compare atom-type numbers and/or atom positions if required.

            if match and (compareAtomTypeNumbers or compareAtomPositions):
                # First check whether the structures have the same number of atoms.

                match = compareStructure.GetAtomCount() == atomCount;

                # Compare atom-type numbers if required.

                if compareAtomTypeNumbers:
                    atomTypeNumbers1 = structure.GetAtomTypeNumbersNumPy(copy = False);
                    atomTypeNumbers2 = compareStructure.GetAtomTypeNumbersNumPy(copy = False);

                    match = (atomTypeNumbers1 == atomTypeNumbers2).all();

                # Compare atom positions if required.

                if match and compareAtomPositions:
                    if transformedPositions is None:
                        if symmetryExpansion != 'none':
                            # If performing symmetry expansions, generate symmetry-transformed positions for the new structure.
                            # If using the 'full' setting, we reduce the expansion to unique structures, mainly to keep the memory taken up by the symmetry-expansion cache to a minimum.

                            transformedPositions = _GenerateSymmetryExpandedPositions(
                                structure, parentSymmetryOperations, tolerance, reduceStructures = symmetryExpansion == 'full'
                                );
                        else:
                            # If not, just compare positions directly.

                            transformedPositions = np.array(
                                [structure.GetAtomPositionsNumPy(copy = False)]
                                );

                    match = _CompareAtomPositions(
                        compareStructure.GetAtomPositionsNumPy(copy = False), transformedPositions, tolerance,
                        compareAtomIndexRanges = compareAtomIndexRanges if compareAtomIndexRanges != None else [(0, atomCount)]
                        );

                    if not match and symmetryExpansion == 'full':
                        # If performing 'full' symmetry expansion, compare tranformations of the structure in the set to tthe ransformations of the new structure.

                        compareTransformedPositions = symmetryExpansionsList[i];

                        if compareTransformedPositions is None:
                            # If the cached symmetry expansions of the reference structure have not been initialised, do so.

                            compareTransformedPositions = _GenerateSymmetryExpandedPositions(
                                compareStructure, parentSymmetryOperations, tolerance, reduceStructures = True
                                );

                            symmetryExpansionsList[i] = compareTransformedPositions;

                        # Compare the positions from each transformation of the reference structure to the transformations of the new structure.

                        for comparePositions in compareTransformedPositions:
                            match = _CompareAtomPositions(
                                comparePositions, transformedPositions, tolerance,
                                compareAtomIndexRanges = compareAtomIndexRanges if compareAtomIndexRanges != None else [(0, atomCount)]
                                );

                            if match:
                                break;

            if match:
                # If the new structure matches, record the index of the match and break.

                matchIndex = i;
                break;

        if matchIndex != None:
            # If a matching structure was found in the set, return a (spacegroup, structure_number) tuple specifying the key.

            return (structure.GetSpacegroup(), matchIndex);
        else:
            return None;

    def _UpdateStructureSet(self, key, structure, degeneracy):
        """
        Update the internal structure set.

        Arguments:
            key -- (spacegroup, structure_number) index of a matching structure in the internal structure set (if found) returned by the _FindStructure() routine.
            structure -- structure to add/update.
            degeneracy -- degeneracy.
        """

        structureSet = self._structureSet;
        symmetryExpansionsCache = self._symmetryExpansionsCache;

        if key is not None:
            # structure is already in the internal structure set -> update its degeneracy.

            spacegroup, index = key;

            _, degeneracyList = structureSet[spacegroup];

            degeneracyList[index] += degeneracy;

            return False;

        else:
            # Add structure to the set.

            spacegroup = structure.GetSpacegroup();

            if spacegroup in structureSet:
                structureList, degeneracyList = structureSet[spacegroup];

                structureList.append(structure);
                degeneracyList.append(degeneracy);

                if symmetryExpansionsCache != None:
                    symmetryExpansionsCache[spacegroup].append(None);
            else:
                structureSet[spacegroup] = ([structure], [degeneracy]);

                if symmetryExpansionsCache != None:
                    symmetryExpansionsCache[spacegroup] = [None];

            return True;

    def _AddStructures(self, structures, degeneracies, isUnion = False, useMP = False, mpNumProcesses = None):
        """
        Merge a list of structures and degeneracies into the set and return the number of structures added.

        Arguments:
            structures -- list of structures to add to the set.
            degeneracies -- list of degeneracies associated with the structures.

        Keyword arguments:
            isUnion -- if True, assume the structures in the list are unique, and only compare them to structures in the initial set while merging.
            useMP -- if True, perform membership testing in parallel using process-based multithreading (requires isUnion = True, default: False).
            mpNumProcesses -- maximum number of worker processes for useMP = True (default: automatically determined).
        """

        structureSet = self._structureSet;

        # Keep track of the number of structures added to the set.

        addCount = None;

        if isUnion:
            # If performing a union, test the structures and add/update as a batch.

            structureKeys = None;

            if useMP and len(structures) > 1:
                # If useMP is set and we have more than one structure, perform the membership testing in parallel.

                numWorkerProcesses = min(
                    mpNumProcesses, len(structures)
                    );

                mappers = [
                    _AddStructures_FindStructureMapper(self)
                        for i in range(0, numWorkerProcesses)
                    ];

                structureKeys = MultiprocessingHelper.QueueMap(
                    structures, mappers, progressBar = False
                    );

            else:
                structureKeys = [
                    self._FindStructure(structure) for structure in structures
                    ];

            # Update the structure set and count the number of structures added.

            addCount = sum(
                1 if self._UpdateStructureSet(key, structure, degeneracy) else 0
                    for structure, degeneracy, key in zip(structures, degeneracies, structureKeys)
                );

        else:
            # Check for and update one structure at a time.

            addCount = 0;

            for structure, degeneracy in zip(structures, degeneracies):
                key = self._FindStructure(structure);

                if self._UpdateStructureSet(key, structure, degeneracy):
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

    def UpdateUnion(self, structureSet, useMP = False, mpNumProcesses = None):
        """
        Perform a union operation with structureSet and return the number of new structures added.

        Arguments:
            structureSet -- structureSet to update with.

        Keyword arguments:
            useMP -- perform the union in parallel using process-based multithreading (default: False).
            mpNumProcesses -- if useMP is set, specify the maximum number of worker processes (default: automatically determined from the CPU core count).

        Notes:
            A union is only valid when both structure sets use the same parameters for identifying unique structures.
        """

        if structureSet == None:
            raise Excpetion("Error: structureSet cannot be None.");

        # Check the union is valid and issue a warning if not.

        isUnion = self.CompareEquivalenceSettings(structureSet);

        if not isUnion:
            warnings.warn("UpdateUnion() is only valid when the supplied structure set is set to use the same comparison parameters as the calling one.", UserWarning);

        return self._AddStructures(
            *structureSet.GetStructureSetFlat(), isUnion = isUnion, useMP = useMP, mpNumProcesses = mpNumProcesses
            );

    def Remove(self, removeIndices):
        """
        Remove structures from the set.

        Arguments:
            removeIndices -- a dictionary with spacegroups as keys and lists of indices to remove as values specifying which structure(s) to remove.

        Notes:
            The key format matches that of the internal structure set dictionary; the keys and indices can be obtained by inspecting the dictionary using the GetStructureSet() method.
            Negative indices are not supported and will raise an IndexError.
        """

        if removeIndices == None:
            raise Exception("Error: structureIndices cannot be None.");

        structureSet = self._structureSet;
        symmetryExpansionsCache = self._symmetryExpansionsCache;

        for spacegroup, indices in removeIndices.items():
            # Validate spacegroup key and structure indices.

            if spacegroup not in structureSet:
                raise KeyError("Error: One or more spacegroup keys was not found in the internal structure set.");

            structures, degeneracies = structureSet[spacegroup];

            for index in indices:
                if index < 0 or index >= len(structures):
                    raise IndexError("Error: One or more indices marked for removal are out of bounds.");

            # Build a new list of structures excluding those specified by the indices.

            newStructures = [
                structure for i, structure in enumerate(structures)
                    if i not in indices
                ];

            # If the new list of structures is not empty, build a new list of degeneracies and update the structure set.
            # If it is, remove the spacegroup key from the set.

            if len(newStructures) != 0:
                newDegeneracies = [
                    degeneracy for i, degeneracy in enumerate(degeneracies)
                        if i not in indices
                    ];

                structureSet[spacegroup] = (newStructures, newDegeneracies);
            else:
                del structureSet[spacegroup];

            # If a cache of symmetry expansions is being kept, update it.

            if symmetryExpansionsCache != None:
                symmetryExpansionsList = symmetryExpansionsCache[i];

                symmetryExpansionsCache[spacegroup] = [
                    expansion for i, expansion in enumerate(symmetryExpansionsList)
                        if i not in indices
                    ];

    def ClearSymmetryExpansionsCache(self):
        """
        Clear the internal symmetry-expansion cache by dereferencing all currently-held sets of symmetry-expanded positions.

        Notes:
            Whether expansions are being cached depends on the symmetryExpansion setting supplied during construction.
            It may be useful to call this method once a StructureSet has been finalised, in order to free any memory being taken up by the cache.
        """

        symmetryExpansionsCache = self._symmetryExpansionsCache;

        if symmetryExpansionsCache != None:
            for spacegroup, symmetryExpansionsList in symmetryExpansionsCache.items():
                symmetryExpansionsCache[spacegroup] = [None] * len(symmetryExpansionsList);

    def CompareEquivalenceSettings(self, structureSet):
        """ Compare the parameters used to determine structural equality with those of structureSet. """

        # Compare tolerances.

        tolerance = self._tolerance;

        # To compare tolerances, we set a threshold an order of magnitude tighter than the minimum tolerance.

        equivalent = math.fabs(tolerance - structureSet._tolerance) < (min(tolerance, structureSet._tolerance) / 10.0);

        if equivalent:
            # Compare symmetry-expansion setting.

            equivalent = self._symmetryExpansion == structureSet._symmetryExpansion;

        if equivalent:
            # Compare symmetry operations.

            if self._parentSymmetryOperations != None:
                if structureSet._parentSymmetryOperations != None:
                    if len(self._parentSymmetryOperations) != len(structureSet._parentSymmetryOperations):
                        equivalent = False;
                    else:
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

        # Compare settings for comparing lattice vectors, atom-type numbers and atom positions.

        equivalent = equivalent and (self._compareLatticeVectors == structureSet._compareLatticeVectors);
        equivalent = equivalent and (self._compareAtomTypeNumbers == structureSet._compareAtomTypeNumbers);
        equivalent = equivalent and (self._compareAtomPositions == structureSet._compareAtomPositions);

        if equivalent:
            # Compare expected atom counts.

            equivalent = self._expectedAtomCount == structureSet._expectedAtomCount;

            if equivalent:
                # Check comparison index ranges, if set.

                if self._compareAtomIndexRanges != None:
                    if structureSet._compareAtomIndexRanges != None:
                        if len(self._compareAtomIndexRanges) != len(structureSet._compareAtomIndexRanges):
                            equivalent = False;
                        else:
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

    def CloneNew(self, structures = None, degeneracies = None, noInitialMerge = False):
        """
        Return a new StructureSet object with the same comparison settings.
        Optionally, pass the structures, degeneracies and noInitialMerge argument to the constructor.
        """

        return StructureSet(
            structures = structures, degeneracies = degeneracies, noInitialMerge = noInitialMerge,
            tolerance = self._tolerance, symmetryExpansion = self._symmetryExpansion, parentSymmetryOperations = self._parentSymmetryOperations,
            compareLatticeVectors = self._compareLatticeVectors, compareAtomTypeNumbers = self._compareAtomTypeNumbers, compareAtomPositions = self._compareAtomPositions,
            expectedAtomCount = self._expectedAtomCount, compareAtomIndexRanges = self._compareAtomIndexRanges
            );


# ----------------------------------------
# _AddStructures_FindStructureMapper Class
# ----------------------------------------

""" Implementation of the MultiprocessingHelper.MapperBase class used for the process-based multithreading in the StructureSet.UpdateUnion method. """

class _AddStructures_FindStructureMapper(MultiprocessingHelper.MapperBase):
    """ Class constructor. """

    def __init__(self, structureSet):
        self._structureSet = structureSet;

    """ Implementation of the MapperBase.Map() method. """

    def Map(self, structure):
        return self._structureSet._FindStructure(structure);


# ------------------
# Internal Functions
# ------------------

def _CompareAtomPositions(comparePositions, refTransformedPositions, tolerance, compareAtomIndexRanges):
    """
    Return True if the atomic positions in comparePositions match any of the set of reference symmetry-transformed positions in refTransformedPositions to within the supplied tolerance.
    Only the position range(s) specified in compareAtomIndexRanges are checked during the comparison.

    Arguments:
        comparePositions -- an Nx3 NumPy matrix containing a set of atom positions.
        refTransformedPositions -- an MxNx3 NumPy matrix containing a set of reference symmetry-transformed atom positions.
        tolerance -- symmetry tolerance for comparing atom positions.
        compareAtomIndexRanges -- range(s) of atom indices to compare.

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

def _GenerateSymmetryExpandedPositions(structure, symmetryOperations, tolerance, reduceStructures = True):
    """
    Generate M sets of symmetry-transformed positions from structure using the supplied symmetry operations.

    Arguments:
        structure -- a Structure object with the positions to be transformed.
        symmetryOperations -- a list of (rotation, translation) tuples specifying the symmetry transformations to apply.
        tolerance -- symmetry tolerance for adjusting transformed positions to minimise numerical noise.

    Keyword arguments:
        reduceStructures -- if True, reduce the expanded structure set to the unique structures.

    Return value:
        An MxNx3 NumPy matrix containing M transformed sets of N atomic positions specified as NumPy vectors.
        The positions are returned sorted by atom type then position, i.e. they match the internal data layout used in the Structure class.
    """

    expandedStructures = None;

    # If the Cython-optimised routines in the _StructureSet module are available, use those.

    if _Cython:
        expandedStructures = _StructureSet._ExpandStructure(structure, symmetryOperations, tolerance);
    else:
        atomData = structure.GetAtomDataNumPy(copy = False);

        numAtoms = len(atomData);
        numSymOps = len(symmetryOperations);

        expandedStructures = np.zeros(
            (numSymOps, numAtoms), dtype = Structure._AtomDataType
            );

        expandedStructures[:] = atomData;

        transformedPositions = expandedStructures.view(dtype = np.float64).reshape((numSymOps, numAtoms, 4))[:, :, 1:];

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
            expandedStructures[i].sort();

    # If required, reduce the set of expanded structures.

    if reduceStructures:
        expandedStructures = np.unique(expandedStructures, axis = 0);

    # Return a copy of the transformed positions.

    return np.copy(
        expandedStructures.view(dtype = np.float64).reshape((len(expandedStructures), structure.GetAtomCount(), 4))[:, :, 1:]
        );
