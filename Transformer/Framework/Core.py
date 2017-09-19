# Transformer/Framework/Core.py by J. M. Skelton


# -------
# Imports
# -------

import multiprocessing;
import warnings;

from Transformer import Structure;
from Transformer import Utility;

from Transformer import _Utility;

from Transformer import StructureSet;


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
        printProgressUpdate = True
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
            Utility.GroupStructuresBySpacegroup(currentStructures, currentDegeneracies, tolerance = tolerance)
            );

    # If symmetryMerge is set, use the symmetry operations of the parent structure during merging.

    parentSymmetryOperations = None;

    if symmetryMerge:
        parentSymmetryOperations = parentStructure.GetSymmetryOperations(tolerance = tolerance);

    # Keep track of the expected number of permutations expected at each substitution.

    permutationCounts = [1];

    # Get a set of atom indices that will be involved in the substitution.

    substitutionAtoms = set();

    for atomType1, atomType2 in atomicSubstitutions:
        substitutionAtoms.add(
            Structure.AtomTypeToAtomTypeNumber(atomType1)
            );

        substitutionAtoms.add(
            Structure.AtomTypeToAtomTypeNumber(atomType2)
            );

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

        if printProgressUpdate:
            print("AtomicSubstitutions(): Performing substitution {0} ({1} -> {2})".format(i + 1, atomType1, atomType2));
            print("AtomicSubstitutions(): Initial structure set contains {0} structure(s)".format(len(currentStructures)));

        # Generate substituted child structures.

        if progressBar:
            print("");

        (newStructures, newDegeneracies), (numGen, numUnique) = _GenerateSubstitutedStructutes(currentStructures, currentDegeneracies, (atomTypeNumber1, atomTypeNumber2), tolerance, substitutionAtoms, parentSymmetryOperations, progressBar);

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

            spacegroupGroups = Utility.GroupStructuresBySpacegroup(
                newStructures, newDegeneracies, tolerance = tolerance
                );

            # If printProgressUpdate is set, print a summary of the spacegroupGroups.

            if printProgressUpdate:
                Utility.PrintSpacegroupGroupSummary(spacegroupGroups);

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
            compareAtomIndexRanges = [atomIndexRanges[atomTypeNumber] for atomTypeNumber in substitutionAtoms if atomTypeNumber in atomIndexRanges];

            structureSet = StructureSet.StructureSet(
                compareLatticeVectors = False, compareAtomTypeNumbers = False,
                tolerance = tolerance, parentSymmetryOperations = parentSymmetryOperations, compareAtomIndexRanges = None #compareAtomIndexRanges
                );

        # Update structureSet with the new structures and degeneracies.

        addCount = structureSet.Update(newStructures, newDegeneracies);

        # Update the counts.

        numGen += len(newStructures);
        numUnique += addCount;

    # Return the merged substituted structures and associated degeneracies along with the numbers of generated/unique substituted structures.

    return ((structureSet.GetStructures(), structureSet.GetDegeneracies()), (numGen, numUnique));


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
