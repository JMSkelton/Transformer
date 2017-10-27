# Transformer/Framework/AtomicSubstitutions.py


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
        parentStructure, substitutions,
        storeIntermediate = None, tolerance = None,
        symmetryMerge = True,
        filterObj = None,
        printProgressUpdate = True,
        useMP = False, mpNumProcesses = None
        ):

    if storeIntermediate != None:
        # If storeIntermediate is provided, sanity-check the indices.

        for index in storeIntermediate:
            if index > len(substitutions):
                raise Exception("Error: If provided, indices in storeIntermediate must be between 0 and the number of substitutions specified by atomicSubstitutions.");
    else:
        # If not initialise it to an index array containing 0 (initial structure) and 1 .. N (all substitutions).

        storeIntermediate = [0] + [i + 1 for i in range(0, len(substitutions))];

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

    # If a filter is provided via the filterObj parameter, perform initialisation.

    if filterObj != None:
        # Check whether the filter is MP safe and, if not, reset useMP and issue a RuntimeWarning.

        if useMP and not filterObj.IsMPSafe():
            warnings.warn("The provided filter is not safe for use with multiprocessing -> useMP will be reset.", RuntimeWarning);

        filterObj.Initialise(substitutions, tolerance, printProgressUpdate, useMP, mpNumProcesses);

    # Keep track of the expected number of permutations expected at each substitution.

    permutationCounts = [1];

    # Decide whether to display a tqdm progress bar while generating child structures.

    progressBar = printProgressUpdate and _TQDM;

    # Perform substitutions.

    for i, substitution in enumerate(substitutions):
        # Get the atom-type numbers of the atoms to find and replace.

        atomType1, atomType2 = substitution;

        atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atomType1);
        atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atomType2);

        # Sanity check.

        if atomTypeNumber1 == None:
            raise Exception("Error: The atom-type number/symbol of the atom to substitute cannot be set to None.");

        # Get indices of atoms that have been manipulated in substitutions up to and including the current one.

        substitutionAtoms = set();

        for atomType1, atomType2 in substitutions[:i + 1]:
            substitutionAtoms.add(
                Structure.AtomTypeToAtomTypeNumber(atomType1)
                );

            substitutionAtoms.add(
                Structure.AtomTypeToAtomTypeNumber(atomType2)
                );

        if printProgressUpdate:
            print("AtomicSubstitutions(): Performing substitution {0} ({1} -> {2})".format(i + 1, atomType1, atomType2));
            print("AtomicSubstitutions(): Initial structure set contains {0} structure(s)".format(len(currentStructures)));

        # If a filter is being used, call its OnStartSubstitution() method with the current substitution index.

        if filterObj != None:
            filterObj.OnStartSubstitution(i);

        # Generate substituted child structures.

        structureSet, numGen = None, None;

        # If useMP is set, use the parallel generation routine; if not, use the serial one.

        if useMP:
            structureSet, numGen = _GenerateSubstitutedStructutesMP(
                currentStructures, currentDegeneracies, (atomTypeNumber1, atomTypeNumber2),
                tolerance, substitutionAtoms, parentSymmetryOperations, filterObj,
                progressBar, mpNumProcesses
                );
        else:
            structureSet, numGen = _GenerateSubstitutedStructutes(
                currentStructures, currentDegeneracies, (atomTypeNumber1, atomTypeNumber2),
                tolerance, substitutionAtoms, parentSymmetryOperations, filterObj,
                progressBar
                );

        if printProgressUpdate:
            numStructures = structureSet.GetStructureCount();

            if numStructures < numGen:
                print("AtomicSubstitutions(): Substituted set contained {0} structure(s)".format(numGen));
                print("AtomicSubstitutions(): Filtering/merging removed {0} structure(s)".format(numGen - numStructures));

                print("");

        # If printProgressUpdate is set, print a summary of the spacegroupGroups.

        if printProgressUpdate:
            StructureTools.PrintSpacegroupGroupSummary(
                structureSet.GetStructureSet()
                );

        # Store the result if required.

        if i + 1 in storeIntermediate:
            intermediateStructures.append(
                structureSet.GetStructureSet()
                );

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

        currentStructures, currentDegeneracies = structureSet.GetStructureSetFlat();

    # If printProgressUpdate is set, print a final summary.

    if printProgressUpdate:
        print("AtomicSubstitutions(): All substitutions performed.");
        print("")

        # Prepend [None] to substitutions for printing the "zeroth" operation (initial structure).

        _PrintResultSummary(
            [None] + substitutions, intermediateStructures, permutationCounts, storeIntermediate
            );

    # After all substitutions have been performed, currentStructures and currentDegeneracies contain the result of the last substitution, and intermediateStructures contains the intermediate results at each step grouped by spacegroup.

    return ((currentStructures, currentDegeneracies), intermediateStructures);


# -----------------------------------------
# Substituted Structure Generation Routines
# -----------------------------------------

def _GenerateSubstitutedStructutes(
        parentStructures, parentDegeneracies, substitution,
        tolerance, substitutionAtoms, parentSymmetryOperations, filterObj,
        progressBar
        ):

    # Unpack the atom-type numbers of the atoms to find and replace.

    atomTypeNumber1, atomTypeNumber2 = substitution;

    # Keep track of how many structures were generated.

    numGen = 0;

    # Set up a primary iterator.

    iValues = range(0, len(parentStructures));

    # If progressBar is set, the iterator is wrapped in the tqdm function, which, if the tqdm module is available, will display a progress bar.

    if progressBar:
        # Padding before progress bar.

        print("");

        iValues = tqdm.tqdm(iValues);

    # We cannot initialise structure set(s) until we have a first child structure.

    structureSet, structureSetFiltered = None, None;

    for i in iValues:
        # Fetch the parent structure and degeneracies.

        structure, degeneracy = parentStructures[i], parentDegeneracies[i];

        # Generate substituted child structures.

        newStructures, newDegeneracies = structure.GetUniqueAtomicSubstitutions(atomTypeNumber1, atomTypeNumber2, tolerance = tolerance);

        # For book-keeping purposes, the degeneracies obtained from the Structure.GetUniqueAtomicSubstitutions() function need to be multiplied by that of the parent structure.

        newDegeneracies = [
            newDegeneracy * degeneracy for newDegeneracy in newDegeneracies
            ];

        # Update the count.

        numGen += len(newStructures);

        # Initialise structure set(s), if not already done.

        if structureSet == None:
            atomCount = newStructures[0].GetAtomCount();

            # As a performance optimisation, we only compare the positions of atoms involved in the substitution.

            atomIndexRanges = newStructures[0].GetAtomIndexRanges();

            compareAtomIndexRanges = [
                atomIndexRanges[atomTypeNumber] for atomTypeNumber in substitutionAtoms
                    if atomTypeNumber in atomIndexRanges
                ];

            structureSet = StructureSet.StructureSetOpt(
                atomCount, tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = compareAtomIndexRanges
                );

            # If a filter has been supplied and requires the filtered structures, set up a second structure set to store them.

            if filterObj != None and filterObj.RequiresFilteredStructures():
                structureSetFiltered = StructureSet.StructureSetOpt(
                    atomCount, tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = compareAtomIndexRanges
                    );

        # If a filter has been supplied, partition the structures between the two structure sets.
        # If not, merge the new structures into the structure set.

        if filterObj != None:
            # New lists to store the structures passed and rejected by the filter (if required).

            passedStructures, passedDegeneracies = [], [];

            rejectedStructures, rejectedDegeneracies = None, None;

            if structureSetFiltered != None:
                rejectedStructures, rejectedDegeneracies = [], [];

            # Pass each structure through the filter and sort.

            for newStructure, newDegeneracy in zip(newStructures, newDegeneracies):
                if filterObj.TestSubstitutedStructure(newStructure, newDegeneracy):
                    passedStructures.append(newStructure);
                    passedDegeneracies.append(newDegeneracy);

                elif structureSetFiltered != None:
                    rejectedStructures.append(newStructure);
                    rejectedDegeneracies.append(newDegeneracy);

            # Update structure set(s).

            structureSet.Update(passedStructures, passedDegeneracies);

            if structureSetFiltered != None:
                structureSetFiltered.Update(rejectedStructures, rejectedDegeneracies);
        else:
            # Update structureSet with the new structures and degeneracies.

            structureSet.Update(newStructures, newDegeneracies);

    if progressBar:
        # Padding after progress bar.

        print("");

    # If a filter has been supplied, call the SetFilteredStructures() and FinaliseMergedStructureSet() methods.

    if filterObj != None:
        # If required, pass SetFilteredStructures() the list of structures removed by TestSubstitutedStructure() and the associated degeneracies.

        if filterObj.RequiresFilteredStructures():
            filterObj.SetFilteredStructureSet(structureSetFiltered);

        # Pass the structure set to FinaliseMergedStructureSet().

        filterObj.FinaliseMergedStructureSet(structureSet);

    # Return the merged substituted structures, associated degeneracies and the number of generated substituted structures.

    return (structureSet, numGen);

# Polling delay for synchronisation between main and worker processes.

_GenerateSubstitutedStructutesMP_PollDelay = 0.1;

def _GenerateSubstitutedStructutesMP(
        parentStructures, parentDegeneracies, substitution,
        tolerance, substitutionAtoms, parentSymmetryOperations, filterObj,
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

    # If the number of parent structures is < mpNumProcesses, there is no point in creating idle processes.

    numWorkerProcesses = min(mpNumProcesses, len(parentStructures));

    # Initialise worker processes.

    processArgs = (substitution, (tolerance, substitutionAtoms, parentSymmetryOperations, filterObj), (inputQueue, progressCounter, terminateFlag, outputQueue));

    workerProcesses = [
        multiprocessing.Process(target = _GenerateSubstitutedStructutesMP_GenerateProcessMain, args = processArgs)
            for i in range(0, numWorkerProcesses)
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
        # Padding before progress bar.

        print("");

        iValues = tqdm.tqdm(iValues);

    for i in iValues:
        while True:
            numProcessed = progressCounter.Current();

            if numProcessed > i:
                break;

            time.sleep(_GenerateSubstitutedStructutesMP_PollDelay);

    if progressBar:
        # Padding after progress bar.

        print("");

    # Signal the worker processes to return their structure sets.

    terminateFlag.value = 1;

    # Collect the output from each of the worker processes.

    results = [];

    for i in range(0, numWorkerProcesses):
        while True:
            # Try to fetch a result from the output queue; if one is not available, sleep for a delay and try again.

            try:
                structureSetRec, structureSetFilteredRec, numGenRec = outputQueue.get();

                # If a worker thread did not obtain any structures to work on, structureSetRec will be None.

                if structureSetRec != None:
                    results.append(
                        (structureSetRec, structureSetFilteredRec, numGenRec)
                        );

                break;

            except Empty:
                time.sleep(_GenerateSubstitutedStructutesMP_PollDelay);

    # Terminate the worker processes if required.

    for workerProcess in workerProcesses:
        if workerProcess.is_alive():
            workerProcess.terminate();

    # Merge result sets using a parallel "divide and conquer" reduction.

    structureSet, structureSetFiltered, numGen = _GenerateSubstitutedStructutesMP_Reduce(
        results, progressBar
        );

    # If a filter has been supplied, call the SetFilteredStructures() and FinaliseMergedStructureSet() methods.

    if filterObj != None:
        if filterObj.RequiresFilteredStructures():
            filterObj.SetFilteredStructureSet(structureSetFiltered);

        filterObj.FinaliseMergedStructureSet(structureSet);

    # Return the structures and degeneracies from the merged set, along with the number of structures generated.

    return (structureSet, numGen);

def _GenerateSubstitutedStructutesMP_GenerateProcessMain(*args):
    # Unpack arguments.

    substitution, (tolerance, substitutionAtoms, parentSymmetryOperations, filterObj), (inputQueue, progressCounter, terminateFlag, outputQueue) = args;

    atomTypeNumber1, atomTypeNumber2 = substitution;

    # Keep track of the number of structures generated.

    numGen = 0;

    # The first structure is required to obtain the atom index ranges needed to initialise the structure sets.

    structureSet, structureSetFiltered = None, None;

    while True:
        # Try to fetch a structure and degeneracy from the input queue; if one is not available, sleep for a delay and try again.

        try:
            structure, degeneracy = inputQueue.get(block = False);

            # Generate substituted structures and modify degeneracies.

            newStructures, newDegeneracies = structure.GetUniqueAtomicSubstitutions(atomTypeNumber1, atomTypeNumber2, tolerance = tolerance);

            newDegeneracies = [
                newDegeneracy * degeneracy for newDegeneracy in newDegeneracies
                ];

            # Update the counter.

            numGen += len(newStructures);

            # Initialise the structure sets if required.

            if structureSet == None:
                atomCount = newStructures[0].GetAtomCount();

                atomIndexRanges = newStructures[0].GetAtomIndexRanges();

                compareAtomIndexRanges = [
                    atomIndexRanges[atomTypeNumber] for atomTypeNumber in substitutionAtoms
                        if atomTypeNumber in atomIndexRanges
                    ];

                structureSet = StructureSet.StructureSetOpt(
                    atomCount, tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = compareAtomIndexRanges
                    );

                if filterObj != None and filterObj.RequiresFilteredStructures():
                    structureSetFiltered = StructureSet.StructureSetOpt(
                        atomCount, tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = compareAtomIndexRanges
                        );

            # Merge the new structures into the structure set(s), depending on whether or not a filter has been supplied and its requirements.

            if filterObj != None:
                passedStructures, passedDegeneracies = [], [];

                rejectedStructures, rejectedDegeneracies = None, None;

                if structureSetFiltered != None:
                    rejectedStructures, rejectedDegeneracies = [], [];

                for newStructure, newDegeneracy in zip(newStructures, newDegeneracies):
                    if filterObj.TestSubstitutedStructure(newStructure, newDegeneracy):
                        passedStructures.append(newStructure);
                        passedDegeneracies.append(newDegeneracy);

                    elif structureSetFiltered != None:
                        rejectedStructures.append(newStructure);
                        rejectedDegeneracies.append(newDegeneracy);

                structureSet.Update(passedStructures, passedDegeneracies);

                if structureSetFiltered != None:
                    structureSetFiltered.Update(rejectedStructures, rejectedDegeneracies);
            else:
                structureSet.Update(newStructures, newDegeneracies);

            # Increment the progress counter.

            progressCounter.Increment();

        except Empty:
            # If the shared terminateFlag is set by the main process, put the structure set(s) and the counter in the output queue and break to allow the process to terminate.

            if terminateFlag.value == 1:
                outputQueue.put(
                    (structureSet, structureSetFiltered, numGen)
                    );

                break;

            # Sleep for _GenerateSubstitutedStructutesMP_PollDelay before checking the queue again.

            time.sleep(_GenerateSubstitutedStructutesMP_PollDelay);

def _GenerateSubstitutedStructutesMP_Reduce(results, progressBar):
    # If there is only one result in the results list, simply return it.

    if len(results) == 1:
        return results[0];

    # Determine whether we are also reducing sets of filtered structures.

    _, structureSetFiltered, _ = results[0];

    reducingFilteredStructureSets = structureSetFiltered != None;

    # Setting a tqdm-based progress bar for the reduction would be (a) fiddly, and (b) not particularly informative; if a progress bar is requested, we print a set of status messages instead.

    if progressBar:
        numStructures = sum(
            structureSet.GetStructureCount() for structureSet, _, _ in results
            );

        numStructuresFiltered = None;

        if reducingFilteredStructureSets:
            numStructuresFiltered = sum(
                structureSetFiltered.GetStructureCount() for _, structureSetFiltered, _ in results
                );

        if reducingFilteredStructureSets:
            print("AtomicSubstitutions(): Reducing {0} structure set(s) w/ {1} + {2} structure(s)".format(len(results), numStructures, numStructuresFiltered));
        else:
            print("AtomicSubstitutions(): Reducing {0} structure set(s) w/ {1} structure(s)".format(len(results), numStructures));

    # The reduction is done via a divide-and-conquer process where N // 2 reductions are performed at each step.
    # The final reduction is a serial merge, while prior steps can be done in parallel using a process pool.

    # Initialise a process pool if required.
    # The number of worker processes we need is equal to the size of the first reduction group, and will always be at most half the value of mpNumProcesses supplied to the calling _GenerateSubstitutedStructutesMP() function.

    processPool = None;

    numGroups = len(results) // 2;

    if numGroups > 1:
        processPool = multiprocessing.Pool(processes = numGroups);

    # Perform the reduction.

    while len(results) > 1:
        # Record the initial numbers of structures.

        numStructures1 = sum(
            structureSet.GetStructureCount() for structureSet, _, _ in results
            );

        numStructuresFiltered1 = None;

        if reducingFilteredStructureSets:
            numStructuresFiltered1 = sum(
                structureSetFiltered.GetStructureCount() for _, structureSetFiltered, _ in results
                );

        # Calculate the number of reduction groups.

        numGroups = len(results) // 2;

        if numGroups <= 1:
            # Final serial merge.

            # Close the process pool if required.

            if processPool != None:
                processPool.close();

            # Merge the structure sets in results into one, and total the number of generated structures.

            structureSet, structureSetFiltered, numGen = results[0];

            for structureSet2, structureSetFiltered2, numGen2 in results[1:]:
                structureSet.UpdateUnion(structureSet2);

                # To avoid keeping unused structures in memory, delete the reference to the old structure set.

                del structureSet2;

                if structureSetFiltered != None:
                    structureSetFiltered.UpdateUnion(structureSetFiltered2);

                    del structureSetFiltered2;

                numGen += numGen2;

            results = [(structureSet, structureSetFiltered, numGen)];
        else:
            # Group the structure sets into numGroups pairs.

            mapGroups = [(results[i], results[i + 1]) for i in range(0, 2 * numGroups, 2)];

            # Reduce each pair using Pool.map().

            resultsNew = processPool.map(
                _GenerateSubstitutedStructutesMP_Reduce_MapFunction, mapGroups
                );

            # Update the results list.

            results = resultsNew + results[2 * numGroups:];

        # If required, print a status message.

        if progressBar:
            numStructures2 = sum(
                structureSet.GetStructureCount() for structureSet, _, _ in results
                );

            numStructuresFiltered2 = None;

            if reducingFilteredStructureSets:
                numStructuresFiltered2 = sum(
                    structureSetFiltered.GetStructureCount() for _, structureSetFiltered, _ in results
                    );

            if numStructures2 < numStructures1 or (reducingFilteredStructureSets and numStructuresFiltered2 < numStructuresFiltered):
                if reducingFilteredStructureSets:
                    print("AtomicSubstitutions(): Reduced {0} + {1} -> {2} + {3} structure(s)".format(numStructures1, numStructuresFiltered1, numStructures2, numStructuresFiltered2));
                else:
                    print("AtomicSubstitutions(): Reduced {0} -> {1} structure(s)".format(numStructures1, numStructures2));

    # Once the merging is complete, return the merged structure set(s) and count.

    if progressBar:
        print("");

    return results[0];

def _GenerateSubstitutedStructutesMP_Reduce_MapFunction(args):
    # Unpack arguments.

    (structureSet1, structureSetFiltered1, numGen1), (structureSet2, structureSetFiltered2, numGen2) = args;

    # Merge the second structure set into the first and delete the former.

    structureSet1.UpdateUnion(structureSet2);

    del structureSet2;

    # If a filter is being used, merge the filtered structure sets.

    if structureSetFiltered1 != None:
        structureSetFiltered1.UpdateUnion(structureSetFiltered2);

        del structureSetFiltered2;

    # Return the updated first structure set(s) along with updated count.

    return (structureSet1, structureSetFiltered1, numGen1 + numGen2);


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
