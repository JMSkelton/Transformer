# Transformer/StructureTools.py by J. M. Skelton


# -------
# Imports
# -------

import array;
import multiprocessing;
import time;
import warnings;

import numpy as np;

try:
    # Python 2.x.

    from Queue import Empty;
except ImportError:
    # Python >= 3.x.

    from queue import Empty;

from Transformer import Constants;

from Transformer.Structure import Structure;

# Try to import the tqdm function to display a progress bar in the _MergeStructureSet*() routines.

ImportedTQDM = False;

try:
    from tqdm import tqdm;

    ImportedTQDM = True;
except ImportError:
    warnings.warn("The tqdm module could not be imported -> displaying progress bars in MergeStructureSet() will be disabled.", RuntimeWarning);

    pass;

# Try to import the Cython-optimised _StructureTools module.
# This module provides core inner routines that replace the Python/NumPy implementations in this module.

CythonImports = False;

try:
    import pyximport; pyximport.install(setup_args = { 'include_dirs' : np.get_include() });

    from Transformer import _StructureTools;

    CythonImports = True;
except ImportError:
    warnings.warn("Optimised merging functions require the pyximport module from the Cython package.", RuntimeWarning);

    pass;


# ----------------------
# Public Merging Routine
# ----------------------

def MergeStructureSet(
        structureSet, degeneracies = None, parentSymmetryOperations = None,
        tolerance = None, compareLatticeVectors = None, compareAtomTypes = None,
        progressBar = False,
        useMP = False, mpNumProcesses = None
        ):

    # If the user does not set a tolerance, set it to the default value used by the Structure class.

    if tolerance == None:
        tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

    numStructures = len(structureSet);

    # If degeneracies are not supplied, the variable is initialised to a list of ones.
    # In this case, degeneracies will store a count of the number of occurrences of each structure.

    if degeneracies == None:
        degeneracies = [1] * numStructures;

    # If useMP is set, dispatch to the parallel merging routine; if not, dispatch to the serial one.

    if useMP:
        return _MergeStructureSetMP(
            structureSet, degeneracies, parentSymmetryOperations,
            tolerance, compareLatticeVectors, compareAtomTypes,
            progressBar,
            mpNumProcesses
            );
    else:
        return _MergeStructureSet(
            structureSet, degeneracies, parentSymmetryOperations,
            tolerance, compareLatticeVectors, compareAtomTypes,
            progressBar
            );


# -----------------------------
# Serial Merging Implementation
# -----------------------------

def _MergeStructureSet(structureSet, degeneracies, parentSymmetryOperations, tolerance, compareLatticeVectors, compareAtomTypes, progressBar):
    numStructures = len(structureSet);

    # Array to keep track of which structures have been marked for removal.

    removeIndices = array.array('b', [0] * numStructures);

    # Set up a primary iterator.
    # If progressBar is set and the tqdm module is available, wrap it in the tqdm() function to display a progress bar.

    iValues = range(0, numStructures);

    if ImportedTQDM and progressBar:
        iValues = tqdm(iValues);

    # Loop over reference structures.

    for i in iValues:
        if removeIndices[i] != 1:
            # Mark equivalent structures.

            equivalentStructureIndices = _MergeStructureSet_MarkDuplicateStructures(structureSet, i, parentSymmetryOperations, tolerance, compareLatticeVectors, compareAtomTypes, removeIndices);

            for j in equivalentStructureIndices:
                # Merge degeneracies and update removeIndices.

                degeneracies[i] += degeneracies[j];

                removeIndices[j] = 1;

    # Discard duplicate structures and redundant degeneracies.

    structureSet = [
        structure for i, structure in enumerate(structureSet)
            if removeIndices[i] == 0
        ];

    degeneracies = [
        degeneracy for i, degeneracy in enumerate(degeneracies)
            if removeIndices[i] == 0
        ];

    return (structureSet, degeneracies);

def _MergeStructureSet_MarkDuplicateStructures(structureSet, refIndex, parentSymmetryOperations, tolerance, compareLatticeVectors, compareAtomTypes, removeIndices):
    refStructure = structureSet[refIndex];

    # List of indices of equivalent structures.

    equivalentStructureIndices = [];

    # Symmetry-transformed positions for the reference structure.
    # These are (relatively) expensive to generate and in some cases may not be used, so we use lazy initialisation.

    refTransformedPositions = None;

    # Loop over other structures.

    for i in range(0, len(structureSet)):
        # Ignore structures that have already been marked for removal.

        if i != refIndex and removeIndices[i] != 1:
            compareStructure = structureSet[i];

            remove = True;

            # These tests are ordered from fastest to slowest.

            if compareLatticeVectors:
                remove = refStructure.CompareLatticeVectors(compareStructure);

            if remove and compareAtomTypes:
                remove = refStructure.CompareAtomTypeNumbers(compareStructure);

            if remove:
                # Initialise refTransformedPositions if required.

                if refTransformedPositions is None:
                    if parentSymmetryOperations != None:
                        # If a set of symmetry operations for a parent structure have been supplied, generate a set of symmetry-transformed positions for the reference structure.

                        if CythonImports:
                            refTransformedPositions = _MergeStructureSet_GenerateSymmetryTransformedPositions(
                                refStructure, parentSymmetryOperations
                                );
                        else:
                            refTransformedPositions = _MergeStructureSet_GenerateSymmetryTransformedPositions(
                                    refStructure, parentSymmetryOperations
                                    );
                    else:
                        # If not, compare against the reference-structure positions themselves.

                        refTransformedPositions = np.array(
                            [refStructure.GetAtomPositionsNumPy(copy = False)]
                            );

                if CythonImports:
                    # If available, use the Cython-optimised position comparison routine.

                    remove = _StructureTools._MergeStructureSet_ComparePositions(
                        compareStructure.GetAtomPositionsNumPy(copy = False), refTransformedPositions, tolerance
                        );
                else:
                    compareResult = np.all(
                        np.abs(refTransformedPositions - compareStructure.GetAtomPositionsNumPy(copy = False)) < tolerance, axis = (1, 2)
                        );

                    remove = compareResult.any();

            # If the structure is a duplicate, add its index to the list.

            if remove:
                equivalentStructureIndices.append(i);

    return equivalentStructureIndices;

def _MergeStructureSet_GenerateSymmetryTransformedPositions(structure, symmetryOperations):
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

    # Sort the data blocks.

    for i in range(0, numSymOps):
        transformedStructures[i].sort();

    # Eliminate redundant transformed structures.
    # This is a relatively expensive operation, but appears to be worth it when using the pure-Python merging routines.

    transformedStructures = np.unique(transformedStructures, axis = 0)

    # Return a view to the transformed positions.

    return transformedStructures.view(dtype = np.float64).reshape((len(transformedStructures), structure.GetAtomCount(), 4))[:, :, 1:];


# -------------------------------
# Parallel Merging Implementation
# -------------------------------

_MergeStructureSetMP_PollDelay = 0.1;

def _MergeStructureSetMP(structureSet, degeneracies, parentSymmetryOperations, tolerance, compareLatticeVectors, compareAtomTypes, progressBar, mpNumProcesses):
    numStructures = len(structureSet);

    # If mpProcesses is not set, set it to the number of physical CPUs.

    if mpNumProcesses == None:
        mpNumProcesses = 1;

        # According to the documentation, cpu_count() may raise a NotImplementedError; if this happens, issue a warning.

        try:
            mpNumProcesses = multiprocessing.cpu_count();
        except NotImplementedError:
            warnings.warn("multiprocessing.cpu_count() is not implemented on this platform, so mpNumProcesses will default to 1.", RuntimeWarning);

    # Create the indices of structures marked for removal as a shared-memory array.

    removeIndices = multiprocessing.Array('b', [0] * numStructures);

    # Set up input and output queues for communication between the main thread and worker processes.

    inputQueue, outputQueue = multiprocessing.Queue(), multiprocessing.Queue();

    # Arguments to pass to _MergeStructureSetMP_ProcessMain().

    workerProcessArgs = (
        structureSet, parentSymmetryOperations,
        tolerance, compareLatticeVectors, compareAtomTypes,
        removeIndices, inputQueue, outputQueue
        );

    # Create and start worker processes.

    workerProcesses = [
        multiprocessing.Process(target = _MergeStructureSetMP_ProcessMain, args = workerProcessArgs)
            for i in range(0, mpNumProcesses)
        ];

    for workerProcess in workerProcesses:
        workerProcess.start();

    # The input queue stores indices of structures to compare to the others in the set.
    # These are retrieved by worker processes, which run _MergeStructureSet_MarkDuplicateStructures and place the index and a list of indices of equivalent structures in the output queue.
    # The main thread retrieves these as they become available, updates the shared-memory removeIndices array, and handles merging the degeneracies.
    # The worker threads check the shared-memory array regularly during the loop and abort if the structure they are working on is marked for removal.
    # While this may lead to some duplicated work compared to the serial merging routine, it appears to be made up for by multicore processing.

    # Dispatch workloads to the input queue.

    for i in range(0, numStructures):
        inputQueue.put(i);

    # Set up a primary iterator for retrieving results fro the output queue.
    # If progressBar is set and the tqdm module is available, wrap it in the tqdm() function to display a progress bar.

    iValues = range(0, numStructures);

    if ImportedTQDM and progressBar:
        iValues = tqdm(iValues);

    # Retrieve results from the output queue.

    for i in iValues:
        # Try to retrieve an index and list of equivalent structures.
        # If no results are available, sleep for _MergeStructureSetMP_PollDelay and try again.

        index, removeMarkIndices = None, None;

        while True:
            try:
                index, equivalentStructureIndices = outputQueue.get(block = False);
                break;
            except Empty:
                time.sleep(_MergeStructureSetMP_PollDelay);

        # If the primary index has not been marked for removal and indices were returned, loop through them and merge the degeneracies, and mark the equivalent structures for removal in the shared-memory array.

        if removeIndices[index] == 0 and equivalentStructureIndices != None:
            for j in equivalentStructureIndices:
                degeneracies[index] += degeneracies[j];

                removeIndices[j] = 1;

    # Terminate the worker processes.

    for workerProcess in workerProcesses:
        workerProcess.terminate();

    # Eliminate duplicate structures and redundant degeneracies.

    structureSet = [
        structure for i, structure in enumerate(structureSet)
            if removeIndices[i] == 0
        ];

    degeneracies = [
        degeneracy for i, degeneracy in enumerate(degeneracies)
            if removeIndices[i] == 0
        ];

    return (structureSet, degeneracies);

def _MergeStructureSetMP_ProcessMain(structureSet, parentSymmetryOperations, tolerance, compareLatticeVectors, compareAtomTypes, removeIndices, inputQueue, outputQueue):
    # Loop until the process is terminated.

    while True:
        # Try to pull a primary index to work on from the input queue.
        # If one is not available, sleep for _MergeStructureSetMP_PollDelay and try again.

        refIndex = None;

        while True:
            try:
                refIndex = inputQueue.get(block = False);
                break;
            except Empty:
                time.sleep(_MergeStructureSetMP_PollDelay);

        equivalentStructureIndices = None;

        # Check the reference structure hasn't already been marked for removal.

        if removeIndices[refIndex] != 1:
            # Compare the reference structure to others in the set and mark duplicates.

            equivalentStructureIndices = _MergeStructureSet_MarkDuplicateStructures(
                structureSet, refIndex, parentSymmetryOperations, tolerance, compareLatticeVectors, compareAtomTypes, removeIndices
                );

            # Put the reference index and list of indices of equivalent structures into the output queue.

        outputQueue.put(
            (refIndex, equivalentStructureIndices)
            );


# -----------------
# Utility Functions
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
