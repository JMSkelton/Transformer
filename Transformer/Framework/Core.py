# Transformer/Framework/AtomicSubstitutions.py


# -------
# Imports
# -------

import multiprocessing;
import time;

from Transformer import Structure;
from Transformer import StructureSet;

from Transformer.Framework import Utilities;

from Transformer.Utilities import MultiprocessingHelper;
from Transformer.Utilities import StructureTools;


# ---------------
# Public Routines
# ---------------

def AtomicSubstitutions(
        structureOrStructureSet, substitutions,
        tolerance = None, symmetryExpansion = 'full', parentSymmetryOperations = None,
        filterObj = None,
        printProgressUpdate = True,
        useMP = False, mpNumProcesses = None
        ):

    # Construct a generator to perform the substitutions.

    generator = _AtomicSubstitutionsIter(
        structureOrStructureSet, substitutions,
        tolerance = tolerance, symmetryExpansion = symmetryExpansion, parentSymmetryOperations = parentSymmetryOperations,
        filterObj = filterObj,
        printProgressUpdate = printProgressUpdate,
        useMP = useMP, mpNumProcesses = mpNumProcesses
        );

    # Keep track of the substitutions performed.

    substitutions = [];

    # Store intermediate structure sets and the substitution numbers at which they have been captured.

    storeIntermediate, intermediateStructureSets = [], [];

    # Store the number of permutations at each round of substitution.

    permutationCounts = [];

    # Run the generator to collect the results.

    for i, (substitution, structureSet, permutationCount) in enumerate(generator):
        substitutions.append(substitution);

        if structureSet != None:
            storeIntermediate.append(i);
            intermediateStructureSets.append(structureSet);

        permutationCounts.append(permutationCount);

    # If printProgressUpdate is set, print a final summary.

    if printProgressUpdate:
        print("AtomicSubstitutions: All substitutions performed.");
        print("")

        # Prepend [None] to substitutions for printing the "zeroth" operation (initial structure).

        _PrintResultSummary(
            substitutions, intermediateStructureSets, permutationCounts, storeIntermediate
            );

    # Return the list of captured intermediate structures.

    return intermediateStructureSets;

def AtomicSubstitutionsIter(
    structureOrStructureSet, substitutions,
    tolerance = None, symmetryExpansion = 'full', parentSymmetryOperations = None,
    filterObj = None,
    printProgressUpdate = True,
    useMP = False, mpNumProcesses = None
    ):

    # Construct a generator to perform the substitutions.

    generator = _AtomicSubstitutionsIter(
        structureOrStructureSet, substitutions,
        tolerance = tolerance, symmetryExpansion = symmetryExpansion, parentSymmetryOperations = parentSymmetryOperations,
        filterObj = filterObj,
        printProgressUpdate = printProgressUpdate,
        useMP = useMP, mpNumProcesses = mpNumProcesses
        );

    # Run the generator and yield captured structure sets.

    for _, structureSet, _ in generator:
        if structureSet != None:
            yield structureSet;


# -------------------------------------
# _StructureSetAccumulator Helper Class
# -------------------------------------

class _StructureSetAccumulator(MultiprocessingHelper.AccumulatorBase):
    # -----------
    # Constructor
    # -----------

    def __init__(self, substitution, substitutionAtoms, tolerance, symmetryExpansion, parentSymmetryOperations, filterObj):
        # Store parameters.

        self._substitution = substitution;
        self._substitutionAtoms = substitutionAtoms;

        self._tolerance = tolerance;
        self._symmetryExpansion = symmetryExpansion;
        self._parentSymmetryOperations = parentSymmetryOperations;

        self._filterObj = filterObj

        # Structure sets to store the generated and, depending on whether a filter is supplied, filtered structures.

        self._structureSet = None;
        self._structureSetFiltered = None;

        # Keep track of the number of structures generated.

        self._numGen = 0;

    # -------
    # Methods
    # -------

    def Accumulate(self, item):
        # Unpack input item.

        structure, degeneracy = item;

        # Generate substituted structures and new degeneracies.

        atomTypeNumber1, atomTypeNumber2 = self._substitution;
        tolerance = self._tolerance;

        newStructures, newDegeneracies = structure.GetUniqueAtomicSubstitutions(atomTypeNumber1, atomTypeNumber2, tolerance = tolerance);

        newDegeneracies = [
            newDegeneracy * degeneracy for newDegeneracy in newDegeneracies
            ];

        # Update the generated structure count.

        self._numGen += len(newStructures);

        # Load the filter, if using, and structure sets.

        filterObj = self._filterObj;
        structureSet, structureSetFiltered = self._structureSet, self._structureSetFiltered;

        # Initialise the structure sets if required; we need the first set of substituted structures in order to obtain the atom index ranges needed to initialise the structure sets.

        if structureSet == None:
            substitutionAtoms = self._substitutionAtoms;

            symmetryExpansion = self._symmetryExpansion;
            parentSymmetryOperations = self._parentSymmetryOperations;

            atomCount = newStructures[0].GetAtomCount();

            compareAtomIndexRanges = None;

            if substitutionAtoms != None:
                atomIndexRanges = newStructures[0].GetAtomIndexRanges();

                compareAtomIndexRanges = [
                    atomIndexRanges[atomTypeNumber] for atomTypeNumber in substitutionAtoms
                        if atomTypeNumber in atomIndexRanges
                    ];
            else:
                compareAtomIndexRanges = [(0, atomCount)];

            structureSet = StructureSet.StructureSet(
                tolerance = tolerance, symmetryExpansion = symmetryExpansion, parentSymmetryOperations = parentSymmetryOperations,
                compareLatticeVectors = False, compareAtomTypeNumbers = False, compareAtomPositions = True,
                expectedAtomCount = atomCount, compareAtomIndexRanges = compareAtomIndexRanges
                );

            self._structureSet = structureSet;

            if filterObj != None and filterObj.RequiresFilteredStructures():
                structureSetFiltered = StructureSet.StructureSet(
                    tolerance = tolerance, symmetryExpansion = symmetryExpansion, parentSymmetryOperations = parentSymmetryOperations,
                    compareLatticeVectors = False, compareAtomTypeNumbers = False, compareAtomPositions = True,
                    expectedAtomCount = atomCount, compareAtomIndexRanges = compareAtomIndexRanges
                    );

                self._structureSetFiltered = structureSetFiltered;

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

    def Finalise(self):
        # Return the collected and filtered structure sets and the number of structures generated.

        return (self._structureSet, self._structureSetFiltered, self._numGen);


# ------------------------------------
# Main Substitution Generator Function
# ------------------------------------

def _AtomicSubstitutionsIter(
    structureOrStructureSet, substitutions,
    tolerance = None, symmetryExpansion = 'full', parentSymmetryOperations = None,
    filterObj = None,
    printProgressUpdate = True,
    useMP = False, mpNumProcesses = None
    ):

    # Perform initialisation.

    initialStructureSet, initialPermutationCount, substitutions, storeIntermediate, symmetryExpansion, parentSymmetryOperations, useIndexRanges = _AtomicSubstitutions_Initialise(
        structureOrStructureSet, substitutions,
        tolerance, symmetryExpansion, parentSymmetryOperations,
        filterObj, printProgressUpdate, useMP, mpNumProcesses
        );

    # Yield the initial structure set and permutation count.

    yield (None, initialStructureSet, initialPermutationCount);

    # Perform substitutions.

    currentStructureSet, currentPermutationCount = initialStructureSet, initialPermutationCount;

    for substitutionIndex in range(0, len(substitutions)):
        currentStructureSet, currentPermutationCount = _AtomicSubstitutions_PerformSubstitution(
            currentStructureSet, currentPermutationCount,
            substitutions, substitutionIndex, tolerance, symmetryExpansion, parentSymmetryOperations, filterObj, useIndexRanges,
            printProgressUpdate, useMP, mpNumProcesses
            );

        # Yield the current permutation count, plus the current structure set if it is to be captured.

        if substitutionIndex + 1 in storeIntermediate:
            yield (substitutions[substitutionIndex], currentStructureSet, currentPermutationCount);
        else:
            yield (substitutions[substitutionIndex], None, currentPermutationCount);


# ----------------------
# Initialisation Routine
# ----------------------

def _AtomicSubstitutions_Initialise(
    structureOrStructureSet, substitutions,
    tolerance, symmetryExpansion, parentSymmetryOperations,
    filterObj, printProgressUpdate, useMP, mpNumProcesses
    ):

    if structureOrStructureSet == None:
        raise Exception("Error: structureOrStructureSet cannot be None.");

    # Set up an initial structure set.

    structureSet = None;

    if isinstance(structureOrStructureSet, Structure.Structure):
        # If an initial structure is supplied, create a structure set from it.

        structureSet = StructureSet.StructureSet(
            structures = [structureOrStructureSet], degeneracies = [1], symmetryExpansion = 'none'
            );
    elif isinstance(structureOrStructureSet, StructureSet.StructureSet):
        # If a structure set is supplied, use that.

        structureSet = structureOrStructureSet;
    else:
        raise Exception("Error: structureOrStructureSet must be either a Structure, a StructureSet, or a derived class.");

    structuresFlat, degeneraciesFlat = structureSet.GetStructureSetFlat();

    # If symmetryMerge is set, and parentSymmetryOperations is not, set it automatically if sensible to do so.

    if symmetryExpansion != 'none' and parentSymmetryOperations is None:
        if len(structuresFlat) == 1:
            parentSymmetryOperations = structuresFlat[0].GetSymmetryOperations(tolerance = tolerance);
        else:
            raise Exception("Error: When an initial structure set is supplied, if symmetryExpansion is not 'none' (the default is 'full'), a set of parent symmetry operations must be supplied.");

    # If the initial structure set contains more than one structure, we need to check that all the structures have the same atomic composition.

    if structureSet.GetStructureCount() > 1:
        atomCount = structuresFlat[0].GetAtomCount();
        atomTypeNumbers = structuresFlat[0].GetAtomTypeNumbersNumPy(copy = False);

        for structure in structuresFlat[1:]:
            if structure.GetAtomCount() != atomCount:
                raise Exception("Error: All structures in the initial structure set must contain the same number of atoms.");

            if (structure.GetAtomTypeNumbersNumPy(copy = False) != atomTypeNumbers).any():
                raise Exception("Error: All structures in the initial structure set must have the same atomic composition.");

    # Unpack substitutions and set up a list of indices of intermediate results to store.

    if substitutions == None or len(substitutions) == 0:
        raise Exception("Error: substitutions cannot be None and must contain at least one element.");

    substitutionsFlat = [];
    storeIntermediate = [0];

    for substitution in substitutions:
        # Each entry in substitutions may be a tuple or a list of tuples.

        if isinstance(substitution, list):
            substitutionsFlat += substitution;
        else:
            substitutionsFlat.append(substitution);

        storeIntermediate.append(
            len(substitutionsFlat)
            );

    # Check that the sequence of substitutions can be performed.

    typeNumbers, atomCounts = structuresFlat[0].GetAtomTypeNumbersCounts();

    # Keep a running count of the composition.

    runningCount = {
        typeNumber : atomCount
            for typeNumber, atomCount in zip(typeNumbers, atomCounts)
        };

    sequenceValid = True;

    for atomType1, atomType2 in substitutionsFlat:
        # Convert substitution atom types to atom-type numbers.

        atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atomType1);
        atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atomType2);

        # Check the type number to be substituted is in the current composition.

        if atomTypeNumber1 not in runningCount:
            sequenceValid = False;
            break;

        if runningCount[atomTypeNumber1] == 0:
            sequenceValid = False;
            break;

        # Update the running count.

        runningCount[atomTypeNumber1] -= 1;

        if atomTypeNumber2 != None:
            if atomTypeNumber2 in runningCount:
                runningCount[atomTypeNumber2] += 1;
            else:
                runningCount[atomTypeNumber2] = 1;

    if not sequenceValid:
        raise Exception("Error: Supplied sequence of substitutions is not valid for the composition of the supplied structure(s).");

    # Set the initial permutation count to the sum of the degeneracies.

    initialPermutationCount = sum(degeneraciesFlat);

    # One of the performance optimisations used by the StructureSet class is to restrict comparing atom positions to specific ranges of indices.
    # If we start from a single structure, we can work this out based on the sequence of substitutions.
    # If we start from a set of structures, however, we cannot be sure this is valid, so we disable it.

    useIndexRanges = structureSet.GetStructureCount() == 1;

    # Finally, if a filter is provided via the filterObj keyword argument, initialise it.

    if filterObj != None:
        filterObj.Initialise(substitutions, tolerance, printProgressUpdate, useMP, mpNumProcesses);

    # Return the initial structure set and permutation count, along with the parent symmetry operations (if using).

    return (structureSet, initialPermutationCount, substitutionsFlat, storeIntermediate, symmetryExpansion, parentSymmetryOperations, useIndexRanges);


# --------------------
# Substitution Routine
# --------------------

def _AtomicSubstitutions_PerformSubstitution(
    currentStructureSet, currentPermutationCount,
    substitutions, substitutionIndex, tolerance, symmetryExpansion, parentSymmetryOperations, filterObj, useIndexRanges,
    printProgressUpdate, useMP, mpNumProcesses
    ):

    # Get the atom-type numbers of the atoms to find and replace.

    atomType1, atomType2 = substitutions[substitutionIndex];

    atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atomType1);
    atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atomType2);

    # If useIndexRanges is set, generate a set of the type numbers of atoms that have been manipulated in substitutions up to and including the current one.

    substitutionAtoms = None;

    if useIndexRanges:
        substitutionAtoms = set();

        for atomType1, atomType2 in substitutions[:substitutionIndex + 1]:
            substitutionAtoms.add(
                Structure.AtomTypeToAtomTypeNumber(atomType1)
                );

            substitutionAtoms.add(
                Structure.AtomTypeToAtomTypeNumber(atomType2)
                );

    if printProgressUpdate:
        print("AtomicSubstitutions: Performing substitution {0} ({1} -> {2})".format(substitutionIndex + 1, atomType1, atomType2));
        print("AtomicSubstitutions: Initial structure set contains {0} structure(s)".format(currentStructureSet.GetStructureCount()));

        print("");

    # If a filter is being used, call its OnStartSubstitution() method with the current substitution index.

    if filterObj != None:
        filterObj.OnStartSubstitution(substitutionIndex);

    # Generate substituted child structures.

    newStructureSet, numGen = None, None;

    newStructureSet, numGen = _GenerateSubstitutedStructutes(
        currentStructureSet, (atomTypeNumber1, atomTypeNumber2),
        tolerance, substitutionAtoms, symmetryExpansion, parentSymmetryOperations, filterObj,
        printProgressUpdate, useMP, mpNumProcesses
        );

    # Since we are finished adding structures to the new structure set, we can immediately clear the symmetry-expansion cache (if in use) to free memory.

    newStructureSet.ClearSymmetryExpansionsCache();

    if printProgressUpdate:
        numStructures = newStructureSet.GetStructureCount();

        if numStructures < numGen:
            print("AtomicSubstitutions: Substituted set contained {0} structure(s)".format(numGen));
            print("AtomicSubstitutions: Filtering/merging removed {0} structure(s)".format(numGen - numStructures));

            print("");

     # Calculate a new permutation count.

    newPermutationCount = None;

    # All structures should have the same composition -> take the atom-type numbers from a random one in the set.

    structureSet = newStructureSet.GetStructureSet();

    for spacegroup, (structures, degeneracies) in structureSet.items():
        if len(structures) > 0:
            atomTypeNumbers = structures[0].GetAtomTypeNumbersNumPy(copy = False);

            # The sum needs to be explicitly converted to a Python integer to avoid arithmetic overflow during long and/or complex substitution sequences.

            newPermutationCount = currentPermutationCount * (int((atomTypeNumbers == atomTypeNumber1).sum()) + 1);

            break;

    # If printProgressUpdate is set, print a summary of the structure set.

    if printProgressUpdate:
        StructureTools.PrintStructureSetSummary(newStructureSet);

    # Return the new structure set and permutation count.

    return (newStructureSet, newPermutationCount);


# ----------------------------------------
# Substituted Structure Generation Routine
# ----------------------------------------

def _GenerateSubstitutedStructutes(
        initialStructureSet, substitution,
        tolerance, substitutionAtoms, symmetryExpansion, parentSymmetryOperations, filterObj,
        printProgressUpdate, useMP, mpNumProcesses
        ):

    # Convert the initial structure set to flat structure and degeneracy lists.

    parentStructures, parentDegeneracies = initialStructureSet.GetStructureSetFlat();

    # Set up StructureSetAccumulator objects to generate the substituted structures and collect them into structure sets.

    numAccumulators = 1;

    if useMP:
        # If mpNumProcesses is not set, set it using the CPUCount() routine in the MultiprocessingHelper module.

        if mpNumProcesses == None:
            mpNumProcesses = MultiprocessingHelper.CPUCount();

        numAccumulators = min(mpNumProcesses, len(parentStructures));

    # Initialise accumulators.

    accumulators = [
        _StructureSetAccumulator(substitution, substitutionAtoms, tolerance, symmetryExpansion, parentSymmetryOperations, filterObj)
            for i in range(0, numAccumulators)
        ];

    # QueueAccumulate() automatically selects a serial or parallel code path depending on the number of accumulator objects supplied to it.

    results = MultiprocessingHelper.QueueAccumulate(
        [item for item in zip(parentStructures, parentDegeneracies)], accumulators, progressBar = printProgressUpdate
        );

    # If the number of parent structures is comparable to the number of worker processes, workers may return empty structure sets.
    # To avoid making the reduction routine messy, we filter the result set before passing it to _GenerateSubstitutedStructutes_Reduce().

    results = [
        (structureSet, structureSetFiltered, numGen) for structureSet, structureSetFiltered, numGen in results
            if structureSet != None
        ];

    # Merge structure sets in results.

    newStructureSet = Utilities.MergeStructureSets(
        [structureSet for structureSet, _, _ in results], useMP = True, mpNumProcesses = mpNumProcesses, printProgressUpdate = printProgressUpdate, inPlace = True
        );

    # If required, merge filtered structure sets in results.

    newStructureSetFiltered = None;

    _, structureSetFiltered, _ = results[0];

    if structureSetFiltered != None:
        newStructureSetFiltered = Utilities.MergeStructureSets(
            [structureSetFiltered for _, structureSetFiltered, _ in results], useMP = True, mpNumProcesses = mpNumProcesses, printProgressUpdate = printProgressUpdate, inPlace = True
            );

    # Sum numbers of generated structures.

    numGen = sum(
        numGen for _, _, numGen in results
        );

    # If a filter has been supplied, call the SetFilteredStructures() and FinaliseMergedStructureSet() methods.

    if filterObj != None:
        if filterObj.RequiresFilteredStructures():
            filterObj.SetFilteredStructureSet(newStructureSetFiltered);

        filterObj.FinaliseMergedStructureSet(newStructureSet);

    # Return the structures and degeneracies from the merged set, along with the number of structures generated.

    return (newStructureSet, numGen);


# ------------------
# Printing Functions
# ------------------

def _PrintResultSummary(substitutions, intermediateStructureSets, permutationCounts, storeIntermediate):
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

    intermediateStructureSetsPointer = 0;

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

            structureSet = intermediateStructureSets[intermediateStructureSetsPointer];

            # Count the number of structures and sum of degeneracies.

            structureCount = structureSet.GetStructureCount();

            degeneracySum = sum(
                sum(degeneracies) for _, degeneracies in structureSet.GetStructureSet().values()
                );

            dataRowData = dataRowData + [
                "{0:,}".format(structureCount),
                "{0:,}".format(degeneracySum),
                "{0:,}".format(permutationCounts[i]),
                ];

            intermediateStructureSetsPointer = intermediateStructureSetsPointer + 1;
        else:
            # If not, print only the expected number of permutations.

            dataRowData = dataRowData + ["-", "-", "{0:,}".format(permutationCounts[i])];

        print(dataRowFormatCode.format(*dataRowData));

    print("");
