# Transformer/Framework/Core.py by J. M. Skelton


# -------
# Imports
# -------

import multiprocessing;
import time;
import warnings;

try:
    # Python 2.x.

    from Queue import Empty;
except ImportError:
    # Newer Python versions.

    from queue import Empty;

from Transformer import Structure;
from Transformer import StructureSet;

from Transformer.Utilities import MultiprocessingHelper;
from Transformer.Utilities import StructureTools;


# Try to import the tqdm module to display progress bars in the _GenerateChildStructutes*() routines.

_TQDM = False;

try:
    import tqdm;

    _TQDM = True;
except ImportError:
    pass;


# ----------------------------------
# "Core" AtomicSubstitutions Routine
# ----------------------------------

def AtomicSubstitutions(
        parentStructure, atomicSubstitutions,
        storeIntermediate = None, tolerance = None,
        symmetryMerge = True,
        printProgressUpdate = True,
        useMP = False, mpNumProcesses = None
        ):

    if storeIntermediate != None:
        # If storeIntermediate is provided, sanity-check the indices.

        for index in storeIntermediate:
            if index > len(atomicSubstitutions):
                raise Exception("Error: If provided, indices in storeIntermediate must be between 0 and the number of substitutions specified by atomicSubstitutions.");
    else:
        # If not initialise it to an index array containing 0 (initial structure) and 1 .. N (all substitutions).

        storeIntermediate = [0] + [i + 1 for i in range(0, len(atomicSubstitutions))];

    # Variables to keep track of the current set of structures and associated degeneracies.

    currentStructures = [parentStructure];
    currentDegeneracies = [1];

    # Store intermediate structure sets obtained after applying substitutions if required.

    intermediateStructures = [];

    if 0 in storeIntermediate:
        intermediateStructures.append(
            StructureTools.GroupStructuresBySpacegroup(currentStructures, currentDegeneracies, tolerance = tolerance)
            );

    # If symmetryMerge is set, use the symmetry operations of the parent structure during merging.

    parentSymmetryOperations = None;

    if symmetryMerge:
        parentSymmetryOperations = parentStructure.GetSymmetryOperations(tolerance = tolerance);

    # Keep track of the expected number of permutations expected at each substitution.

    permutationCounts = [1];

    # Decide whether to display a tqdm progress bar while generating child structures.

    progressBar = printProgressUpdate and _TQDM;

    # Perform substitutions.

    for i, substitution in enumerate(atomicSubstitutions):
        # Get the atom-type numbers of the atoms to find and replace.

        atomType1, atomType2 = substitution;

        atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atomType1);
        atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atomType2);

        # Sanity check.

        if atomTypeNumber1 == None:
            raise Exception("Error: The atom-type number/symbol of the atom to substitute cannot be set to None.");

        # Get indices of atoms that have been manipulated in substitutions up to and including the current one.

        substitutionAtoms = set();

        for atomType1, atomType2 in atomicSubstitutions[:i + 1]:
            substitutionAtoms.add(
                Structure.AtomTypeToAtomTypeNumber(atomType1)
                );

            substitutionAtoms.add(
                Structure.AtomTypeToAtomTypeNumber(atomType2)
                );

        if printProgressUpdate:
            print("AtomicSubstitutions(): Performing substitution {0} ({1} -> {2})".format(i + 1, atomType1, atomType2));
            print("AtomicSubstitutions(): Initial structure set contains {0} structure(s)".format(len(currentStructures)));

        # Generate substituted child structures.

        if progressBar:
            print("");

        newStructures, newDegeneracies = None, None;
        numGen, numUnique = None, None;

        # If useMP is set, use the parallel generation routine; if not, use the serial one.

        if useMP:
            (newStructures, newDegeneracies), (numGen, numUnique) = _GenerateSubstitutedStructutesMP(
                currentStructures, currentDegeneracies, (atomTypeNumber1, atomTypeNumber2),
                tolerance, substitutionAtoms, parentSymmetryOperations,
                progressBar, mpNumProcesses
                );
        else:
            (newStructures, newDegeneracies), (numGen, numUnique) = _GenerateSubstitutedStructutes(
                currentStructures, currentDegeneracies, (atomTypeNumber1, atomTypeNumber2),
                tolerance, substitutionAtoms, parentSymmetryOperations,
                progressBar
                );

        if progressBar:
            print("");

        if printProgressUpdate:
            print("AtomicSubstitutions(): Substituted set contained {0} structure(s)".format(numGen));
            print("AtomicSubstitutions(): Merging removed {0} structure(s)".format(numGen - numUnique));
            print("");

        # Check whether a progress update has been requested or we need to store this set of intermediate results.
        # If neither, we may be able to avoid sorting the structures by spacegroup, and hence some spglib calls.

        if printProgressUpdate or i + 1 in storeIntermediate:
            # Group the new structures by spacegroup.

            spacegroupGroups = StructureTools.GroupStructuresBySpacegroup(
                newStructures, newDegeneracies, tolerance = tolerance
                );

            # If printProgressUpdate is set, print a summary of the spacegroupGroups.

            if printProgressUpdate:
                StructureTools.PrintSpacegroupGroupSummary(spacegroupGroups);

            # Store the result if required.

            if i + 1 in storeIntermediate:
                intermediateStructures.append(spacegroupGroups);

        # Update the permutation count; all structures should have the same composition -> take the atom-type numbers from the first one.

        atomTypes = currentStructures[0].GetAtomTypeNumbers();

        atomCount = 0;

        for atomType in atomTypes:
            if atomType == atomTypeNumber1:
                atomCount = atomCount + 1;

        permutationCounts.append(
            permutationCounts[-1] * atomCount
            );

        # Update the current structure set/degeneracies.

        currentStructures, currentDegeneracies = newStructures, newDegeneracies;

    # If printProgressUpdate is set, print a final summary.

    if printProgressUpdate:
        print("AtomicSubstitutions(): All substitutions performed.");
        print("")

        # Prepend [None] to atomicSubstitutions for printing the "zeroth" operation (initial structure).

        _PrintResultSummary(
            [None] + atomicSubstitutions, intermediateStructures, permutationCounts, storeIntermediate
            );

    # After all substitutions have been performed, currentStructures and currentDegeneracies contain the result of the last substitution, and intermediateStructures contains the intermediate results at each step grouped by spacegroup.

    return ((currentStructures, currentDegeneracies), intermediateStructures);


# -----------------------------------------
# Substituted Structure Generation Routines
# -----------------------------------------

def _GenerateSubstitutedStructutes(
        parentStructures, parentDegeneracies, substitution,
        tolerance, substitutionAtoms, parentSymmetryOperations,
        progressBar
        ):

    # Unpack the atom-type numbers of the atoms to find and replace.

    atomTypeNumber1, atomTypeNumber2 = substitution;

    # Keep track of how many structures were generated and added to the merged structure set.

    numGen, numUnique = 0, 0;

    # Set up a primary iterator.
    # If progressBar is set, the iterator is wrapped in the tqdm function, which, if the tqdm module is available, will display a progress bar.

    iValues = range(0, len(parentStructures));

    if progressBar:
        iValues = tqdm.tqdm(iValues);

    # We cannot initialise structureSet until we have a first child structure.

    structureSet = None;

    for i in iValues:
        # Fetch the parent structure and degeneracies.

        structure, degeneracy = parentStructures[i], parentDegeneracies[i];

        # Generate substituted child structures.

        newStructures, newDegeneracies = structure.GetUniqueAtomicSubstitutions(atomTypeNumber1, atomTypeNumber2, tolerance = tolerance);

        # For book-keeping purposes, the degeneracies obtained from the Structure.GetUniqueAtomicSubstitutions() function need to be multiplied by that of the parent structure.

        newDegeneracies = [
            newDegeneracy * degeneracy for newDegeneracy in newDegeneracies
            ];

        # Initialise structureSet, if not already done.

        if structureSet == None:
            # As a performance optimisation, we only compare the positions of atoms involved in the substitution.

            atomIndexRanges = newStructures[0].GetAtomIndexRanges();

            compareAtomIndexRanges = [
                atomIndexRanges[atomTypeNumber] for atomTypeNumber in substitutionAtoms
                    if atomTypeNumber in atomIndexRanges
                ];

            structureSet = StructureSet.StructureSet(
                compareLatticeVectors = False, compareAtomTypeNumbers = False, compareAtomPositions = True,
                tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = compareAtomIndexRanges
                );

        # Update structureSet with the new structures and degeneracies.

        addCount = structureSet.Update(newStructures, newDegeneracies);

        # Update the counts.

        numGen += len(newStructures);
        numUnique += addCount;

    # Return the merged substituted structures and associated degeneracies along with the numbers of generated/unique substituted structures.

    return ((structureSet.GetStructures(), structureSet.GetDegeneracies()), (numGen, numUnique));

# Polling delay for synchronisation between main and worker processes.

_GenerateSubstitutedStructutesMP_PollDelay = 0.1;

def _GenerateSubstitutedStructutesMP(
        parentStructures, parentDegeneracies, substitution,
        tolerance, substitutionAtoms, parentSymmetryOperations,
        progressBar, mpNumProcesses
        ):

    # Queue to pass parent structures to worker processes.

    inputQueue = multiprocessing.Queue();

    # Synchronised counter for worker processes to report progress.

    progressCounter = MultiprocessingHelper.Counter();

    # Shared-memory flag to signal worker processes to terminate.

    terminateFlag = multiprocessing.Value('B', 0);

    # Queue for worker processes to return merged structure sets.

    outputQueue = multiprocessing.Queue();

    # If mpNumProcesses is not set, set it using the CPUCount() routine in the MultiprocessingHelper module.

    if mpNumProcesses == None:
        mpNumProcesses = MultiprocessingHelper.CPUCount();

    # Initialise worker processes.

    processArgs = (substitution, (tolerance, substitutionAtoms, parentSymmetryOperations), (inputQueue, progressCounter, terminateFlag, outputQueue));

    workerProcesses = [
        multiprocessing.Process(target = _GenerateSubstitutedStructutesMP_ProcessMain, args = (processArgs, ))
            for i in range(0, mpNumProcesses)
        ];

    # Start the worker processes.

    for workerProcess in workerProcesses:
        workerProcess.start();

    # Put pairs of structures and degeneracies to the input queue.

    for item in zip(parentStructures, parentDegeneracies):
        inputQueue.put(item);

    # Poll the progress counter until the workers are finished.
    # This is done inside a for loop to allow a tqdm progress bar to be displayed.

    iValues = range(0, len(parentStructures));

    if progressBar:
        iValues = tqdm.tqdm(iValues);

    for i in iValues:
        while True:
            numProcessed = progressCounter.Current();

            if numProcessed > i:
                break;

            time.sleep(_GenerateSubstitutedStructutesMP_PollDelay);

    # Signal the worker processes to return their structure sets.

    terminateFlag.value = 1;

    # Collect and merge the output from the worker processes.

    structureSet = None;
    numGen, numUnique = None, None;

    # Allow a progress bar to be displayed during merging.

    iValues = range(0, len(workerProcesses));

    if progressBar:
        iValues = tqdm.tqdm(iValues);

    for i in iValues:
        while True:
            # Try to fetch a result from the output queue; if one is not available, sleep for a delay and try again.

            try:
                structureSetRec, (numGenRec, numUniqueRec) = outputQueue.get();

                # If a worker thread did not obtain any structures to work on, structureSetRec will be None.

                if structureSetRec != None:
                    if structureSet == None:
                        # First set of data -> set structureSet and numGen/numUnique.

                        structureSet = structureSetRec;
                        numGen, numUnique = numGenRec, numUniqueRec;
                    else:
                        # Additional set of data -> merge into the first set.

                        addCount = structureSet.UpdateUnion(structureSetRec);

                        numGen += numGenRec;

                        # numUnique should be updated after merging into the first structure set.

                        numUnique += addCount;

                break;

            except Empty:
                time.sleep(_GenerateSubstitutedStructutesMP_PollDelay);

    # Terminate the worker processes if required.

    for workerProcess in workerProcesses:
        if workerProcess.is_alive():
            workerProcess.terminate();

    # Return the structures and degeneracies from the merged set, along with the counts.

    return ((structureSet.GetStructures(), structureSet.GetDegeneracies()), (numGen, numUnique));

def _GenerateSubstitutedStructutesMP_ProcessMain(args):
    # Unpack arguments.

    substitution, (tolerance, substitutionAtoms, parentSymmetryOperations), (inputQueue, progressCounter, terminateFlag, outputQueue) = args;

    atomTypeNumber1, atomTypeNumber2 = substitution;

    # Keep track of the number of structures generated and added to the merged structure set.

    numGen, numUnique = 0, 0;

    # The first structure is required to obtain the atom index ranges needed to initialise structureSet.

    structureSet = None;

    while True:
        # Try to fetch a structure and degeneracy from the input queue; if one is not available, sleep for a delay and try again.

        try:
            structure, degeneracy = inputQueue.get(block = False);

            # Initialise structureSet if required.

            if structureSet == None:
                atomIndexRanges = structure.GetAtomIndexRanges();

                compareAtomIndexRanges = [
                    atomIndexRanges[atomTypeNumber] for atomTypeNumber in substitutionAtoms
                        if atomTypeNumber in atomIndexRanges
                    ];

                structureSet = StructureSet.StructureSet(
                    compareLatticeVectors = False, compareAtomTypeNumbers = False, compareAtomPositions = True,
                    tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = compareAtomIndexRanges
                    );

            # Generate substituted structures, modify degeneracies, and add to the structure set.

            newStructures, newDegeneracies = structure.GetUniqueAtomicSubstitutions(atomTypeNumber1, atomTypeNumber2, tolerance = tolerance);

            newDegeneracies = [
                newDegeneracy * degeneracy for newDegeneracy in newDegeneracies
                ];

            addCount = structureSet.Update(newStructures, newDegeneracies);

            # Update the counters.

            numGen += len(newStructures);
            numUnique += addCount;

            # Increment the progress counter.

            progressCounter.Increment();

        except Empty:
            # If the shared terminateFlag is set by the main process, put the structure set and the counters in the output queue and break to allow the process to terminate.

            if terminateFlag.value == 1:
                outputQueue.put(
                    (structureSet, (numGen, numUnique))
                    );

                break;

            # Sleep for _GenerateSubstitutedStructutesMP_PollDelay before checking the queue again.

            time.sleep(_GenerateSubstitutedStructutesMP_PollDelay);


# ------------------
# Printing Functions
# ------------------

def _PrintResultSummary(substitutions, intermediateStructures, permutationCounts, storeIntermediate):
    # Find the largest permutation count (= maximum integer value to be printed), and get the length of the text fields.

    fieldLength = max(
        len("{0:,}".format(max(permutationCounts))), 15
        );

    # Format and print header row.

    headerRowFormatCode = "{{0: ^16}} | {{1: ^{0}}} | {{2: ^{0}}} | {{3: ^{0}}}".format(fieldLength);

    headerRow = headerRowFormatCode.format("Substitution", "# Structures", "sum(Degeneracy)", "Permutations");

    print(headerRow);
    print('-' * len(headerRow));

    # Generate and print table rows.

    # Depending on the indices in storeIntermediate, the data rows will contain different items of data, making generating the formatted rows somewhat complex.
    # This is the main reason for separating out the formatting code into a separate, private function.

    dataRowFormatCode = "{{0: <3}} {{1: <4}} {{2: <2}} {{3: <4}} | {{4: >{0}}} | {{5: >{0}}} | {{6: >{0}}}".format(fieldLength);

    intermediateStructuresPointer = 0;

    for i, substitution in enumerate(substitutions):
        # The first element of substitutions will be None (input structure; no substitution).
        # The remainder will contain the target/substituted elements in each round of substitution.

        dataRowData = [i];

        if substitution == None:
            dataRowData = dataRowData + ["None", "", ""];
        else:
            atomType1, atomType2 = substitution;

            dataRowData = dataRowData + [atomType1, "->", atomType2 if atomType2 != None else "None"];

        # intermediateStructures will contain a set of structures, grouped by spacegroup, for each index in storeIntermediate.

        if i in storeIntermediate:
            # If intermediate structures following the current substitution were captured, display the number of unique structures along the sum of the degeneracies to compare to the expected number of permutations.

            spacegroupGroups = intermediateStructures[intermediateStructuresPointer];

            # Sum up the number of structures in each group.

            structureCount = sum(
                len(structures) for structures, _ in spacegroupGroups.values()
                );

            # Sum the degeneracies of each structure.

            degeneracySum = sum(
                sum(degeneracies) for _, degeneracies in spacegroupGroups.values()
                );

            dataRowData = dataRowData + [
                "{0:,}".format(structureCount),
                "{0:,}".format(degeneracySum),
                "{0:,}".format(permutationCounts[i]),
                ];

            intermediateStructuresPointer = intermediateStructuresPointer + 1;
        else:
            # If not, print only the expected number of permutations.

            dataRowData = dataRowData + ["-", "-", "{0:,}".format(permutationCounts[i])];

        print(dataRowFormatCode.format(*dataRowData));

    print("");
