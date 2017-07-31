# Transformer/ConvenienceFunctions.py by J. M. Skelton


# -------
# Imports
# -------

import os;
import tarfile;
import warnings;

from Transformer import Constants;
from Transformer import IO;
from Transformer import Structure;
from Transformer import StructureTools;


# ---------------------
# Convenience Functions
# ---------------------

def AntisiteDefects(structure, atom1, atom2, numDefects = None, printProgressUpdate = True, tolerance = None):
    # Convert atom1 and atom2 to atom-type numbers.

    atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atom1);
    atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atom2);

    # Count the number of both atoms in the structure.

    atomTypeNumbers = structure.GetAtomTypeNumbers();

    # Sanity check.

    if atomTypeNumber1 == atomTypeNumber2:
        raise Exception("Error: atom1 and atom2 have the same atom-type number - this is most likely an error.");

    if atomTypeNumber1 not in atomTypeNumbers or atomTypeNumber2 not in atomTypeNumbers:
        raise Exception("Error: One of both of atom1/atom2 were not found in the supplied structure.");

    atomCount1, atomCount2 = 0, 0;

    for atomTypeNumber in atomTypeNumbers:
        if atomTypeNumber == atomTypeNumber1:
            atomCount1 = atomCount1 + 1;
        elif atomTypeNumber == atomTypeNumber2:
            atomCount2 = atomCount2 + 1;

    # The maximum number of possible antisite defects the minimum of atomCount1 and atomCount2.

    maxNumDefects = min(atomCount1, atomCount2);

    # If numDefects is provided, check it; if not, set it to maxNumDefects.

    if numDefects != None:
        if numDefects > maxNumDefects:
            raise Exception("Error: If provided, numDefects cannot be greater than the number of either of atom1/atom2 in the supplied structure.");
    else:
        numDefects = maxNumDefects;

    # When making the substitutions, to avoid reversing earlier substitutions when creating multiple defects, we swap atom1 and atom2 for the arbitrary placeholder symbols.

    # Select a pair of atom-type numbers as placeholders.

    placeholder1, placeholder2 = None, None;

    trialTypeNumber = -1;

    while True:
        if trialTypeNumber not in atomTypeNumbers:
            placeholder1 = trialTypeNumber;
            break;

        trialTypeNumber = trialTypeNumber - 1;

    trialTypeNumber = trialTypeNumber - 1;

    while True:
        if trialTypeNumber not in atomTypeNumbers:
            placeholder2 = trialTypeNumber;
            break;

        trialTypeNumber = trialTypeNumber - 1;

    # Collect structures with antisite defects by performing a sequence of substitutions where we swap atom1 -> 1 and atom2 -> 2.

    _, antisiteDefects = AtomicSubstitutions(
        structure, [(atom1, placeholder1), (atom2, placeholder2)] * numDefects,
        storeIntermediate = [i for i in range(0, 2 * numDefects + 1, 2)], printProgressUpdate = printProgressUpdate, tolerance = tolerance
        );

    # Modify the substituted structures in the result set to swap the placeholders for 1 -> atom2 and 2 -> atom1.
    # After doing so, regenerate the spacegroup groupings and check for duplicate structures.
    # [I couldn't decide whether the second step was actually needed, so I did it anyway, to be on the safe side...!]

    for i, spacegroupGroups in enumerate(antisiteDefects[1:]):
        if printProgressUpdate:
            print("AntisiteDefects(): Post processing defect set {0}.".format(i + 1));
            print("");

        # Merge the structures and degeneracies in the spacegroup groups into a "flat" lists.

        structuresFlat, degeneraciesFlat = [], [];

        for structures, degeneracies in spacegroupGroups.values():
            for structure in structures:
                structuresFlat.append(
                    structure.GetAtomSwap(placeholder1, atom2).GetAtomSwap(placeholder2, atom1)
                    );

            degeneraciesFlat = degeneraciesFlat + degeneracies;

        numStructures = len(structuresFlat);

        # Perform a final merge.

        structuresFlat, degeneraciesFlat = StructureTools.MergeStructureSet(
            structuresFlat, degeneraciesFlat,
            parentSymmetryOperations = structure.GetSymmetryOperations(), tolerance = tolerance,
            compareLatticeVectors = False, compareAtomTypes = False
            );

        reductionCount = len(structuresFlat) - numStructures;

        if printProgressUpdate and reductionCount > 0:
            print("AntisiteDefects(): Final merge removed {0} structure(s)".format(reductionCount));
            print("");

        # Regroup the structures.

        spacegroupGroups = StructureTools.GroupStructuresBySpacegroup(structuresFlat, degeneraciesFlat);

        if printProgressUpdate:
            PrintSpacegroupGroupSummary(spacegroupGroups);

        # Update the results.

        antisiteDefects[i + 1] = spacegroupGroups;

    # Return the result.

    return antisiteDefects;

def SolidSolution(structure, atom1, atom2, printProgressUpdate = True, useShortcut = False, tolerance = None):
    # Convert atom1 and atom2 to atom-type numbers.

    atomFind = Structure.AtomTypeToAtomTypeNumber(atom1);
    atomReplace = Structure.AtomTypeToAtomTypeNumber(atom2);

    atomTypeNumbers = structure.GetAtomTypeNumbers();

    # Sanity checks.

    if atomFind not in atomTypeNumbers:
        raise Exception("Error: The atom specified by atomTypeNumber1/atomicSymbol1 was not found in the supplied structure.");

    if atomReplace in atomTypeNumbers:
        raise Exception("Error: the atom specified by atomTypeNumber2/atomicSymbol2 was found in the supplied structure - this is most likely an error.");

    # Count the number of atoms to substitute.

    substitutionCount = 0;

    for atomType in atomTypeNumbers:
        if atomType == atomFind:
            substitutionCount = substitutionCount + 1;

    # If useShortcut is set, we can save some time by only enumerating the structures for up to ~50 % substitution, and generating the rest by swapping atoms.

    substitutions = None;

    if useShortcut:
        substitutions = [(atom1, atom2)] * (substitutionCount // 2 + substitutionCount % 2);
    else:
        substitutions = [(atom1, atom2)] * substitutionCount;

    # Build the set of solid solutions by performing successive atomic substitutions.

    _, solidSolutions = AtomicSubstitutions(
        structure, substitutions,
        printProgressUpdate = printProgressUpdate, tolerance = tolerance
        );

    # If useShortcut is set, generate the structures for the remaining substitutions.

    if useShortcut:
        # Select an atom-type number to use as a placeholder.

        placeholder = None;

        trialTypeNumber = -1;

        while True:
            if trialTypeNumber not in atomTypeNumbers and trialTypeNumber != atomReplace:
                placeholder = trialTypeNumber;
                break;

            trialTypeNumber = trialTypeNumber - 1;

        for i in range(len(substitutions) + 1, substitutionCount):
            # If requested, print a progress update.

            if printProgressUpdate:
                print("SolidSolution(): Inverting substitution set {0} -> {1}".format(substitutionCount - i, i));

            newResultSet = { };

            resultSet = solidSolutions[substitutionCount - i];

            for key, (structures, degeneracies) in resultSet.items():
                newStructures = [];

                for structure in structures:
                    # Temporarily swap atomReplace with placeholder.

                    structure = structure.GetAtomSwap(atomReplace, placeholder);

                    # Swap atomFind with atomReplace.

                    structure = structure.GetAtomSwap(atomFind, atomReplace);

                    # Finally, replace placeholder with atomFind.

                    structure = structure.GetAtomSwap(placeholder, atomFind);

                    # Add the new structure to the list.

                    newStructures.append(structure);

                # Add the new structures to the new result set, under the same spacegroup key, with a copy of the degeneracies.

                newResultSet[key] = (
                    newStructures, list(degeneracies)
                    );

            # Extend the list of solid solutions with the new result set.

            solidSolutions.append(newResultSet);

        if printProgressUpdate:
            print("");

        if printProgressUpdate:
            # Work out the number of permutations at each round of substitution.

            permutationCounts = [1];

            for i in range(0, substitutionCount):
                permutationCounts.append(
                    permutationCounts[i] * (substitutionCount - i)
                    );

            # Print a summary of the new results ("repurposing" the _AtomicSubstitutions_PrintResultSummary() function).

            _AtomicSubstitutions_PrintResultSummary(
                [None] + [(atom1, atom2)] * substitutionCount, solidSolutions, permutationCounts, [i for i in range(0, len(solidSolutions))]
                );

    # Return the result.

    return solidSolutions;

def AtomicSubstitutions(structure, atomicSubstitutions, storeIntermediate = None, printProgressUpdate = True, tolerance = None):
    if storeIntermediate != None:
        # If storeIntermediate is provided, sanity-check the indices.

        for index in storeIntermediate:
            if index > len(atomicSubstitutions):
                raise Exception("Error: If provided, indices in storeIntermediate must be between 0 and the number of substitutions specified by atomicSubstitutions.");
    else:
        # If not initialise it to an index array containing 0 (initial structure) and 1 .. N (all substitutions).

        storeIntermediate = [0] + [i + 1 for i in range(0, len(atomicSubstitutions))];

    # Compute the symmetry operations of the parent structure, needed for merging.

    parentSymmetryOperations = structure.GetSymmetryOperations();

    # Variables to keep track of the current set of structures and associated degeneracies.

    currentStructures = [structure];
    currentDegeneracies = [1];

    # Store intermediate structure sets obtained after applying substitutions if required.

    intermediateStructures = [];

    if 0 in storeIntermediate:
        intermediateStructures.append(
            StructureTools.GroupStructuresBySpacegroup(currentStructures, currentDegeneracies)
            );

    # Keep track of the expected number of permutations expected at each substitution.

    permutationCounts = [1];

    # Perform each substitution in sequence.

    for i, substitution in enumerate(atomicSubstitutions):
        # Get the atom-type numbers of the atoms to find and replace.;

        atomType1, atomType2 = substitution;

        atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atomType1);
        atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atomType2);

        # Sanity check.

        if atomTypeNumber1 == None:
            raise Exception("Error: The atom-type number/symbol of the atom to substitute cannot be set to None.");

        if printProgressUpdate:
            print("AtomicSubstitutions(): Performing substitution {0} ({1} -> {2})".format(i + 1, atomType1, atomType2));
            print("AtomicSubstitutions(): Initial structure set contains {0} structure(s)".format(len(currentStructures)));

        newStructures, newDegeneracies = [], [];

        for structure, degeneracy in zip(currentStructures, currentDegeneracies):
            atomTypeNumbers = structure.GetAtomTypeNumbers();

            # Get the indices and site degeneracies of the unique atoms in the structure.

            uniqueAtomIndices, siteDegeneracies = structure.GetUniqueAtomIndices();

            for atomIndex, siteDegeneracy in zip(uniqueAtomIndices, siteDegeneracies):
                if atomTypeNumbers[atomIndex] == atomTypeNumber1:
                    # Clone the structure.

                    newStructure = structure.Clone();

                    # Perform the substitution.

                    newStructure.SetAtom(atomIndex, atomTypeNumber2);

                    # Add the new structure to the list.

                    newStructures.append(newStructure);

                    # For "book keeping", the degeneracy of the substitution is the degeneracy of the original structure multiplied by the degeneracy of the substituted site.

                    newDegeneracies.append(degeneracy * siteDegeneracy);

        if printProgressUpdate:
            print("AtomicSubstitutions(): Substituted structure set contains {0} structure(s)".format(len(newStructures)));

        # Merge the new structure set.

        numStructures = len(newStructures);

        # Disable comparing lattice vectors and atom types to avoid unnecessary comparisons.

        newStructures, newDegeneracies = StructureTools.MergeStructureSet(
            newStructures, newDegeneracies, parentSymmetryOperations = parentSymmetryOperations, tolerance = tolerance,
            compareLatticeVectors = False, compareAtomTypes = False
            );

        reductionCount = numStructures - len(newStructures);

        if printProgressUpdate and reductionCount > 0:
            print("AtomicSubstitutions(): Merging removed {0} structure(s)".format(reductionCount));

        # Check whether a progress update has been requested or we need to store this set of intermediate results.
        # If neither, we can avoid sorting the structures by spacegroup, and hence a potentially-large number of spglib calls.

        if printProgressUpdate or i + 1 in storeIntermediate:
            # Group the new structures by spacegroup.

            spacegroupGroups = StructureTools.GroupStructuresBySpacegroup(newStructures, newDegeneracies);

            # If printProgressUpdate is set, print a summary of the spacegroupGroups.

            if printProgressUpdate:
                print("");

                PrintSpacegroupGroupSummary(spacegroupGroups);

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

        _AtomicSubstitutions_PrintResultSummary([None] + atomicSubstitutions, intermediateStructures, permutationCounts, storeIntermediate);

    # After all substitutions have been performed, currentStructures and currentDegeneracies contain the result of the last substitution, and intermediateStructures contains the intermediate results at each step grouped by spacegroup.

    return ((currentStructures, currentDegeneracies), intermediateStructures);

def _AtomicSubstitutions_PrintResultSummary(substitutions, intermediateStructures, permutationCounts, storeIntermediate):
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


# ----------------------
# Batch Export Functions
# ----------------------

SupportedExportFileFormats = ['vasp', 'aims'];

def ExportAtomicSubstitutionResultSet(resultSet, prefix = None, atomicSymbolLookupTable = None, workingDirectory = "./", fileFormat = 'vasp', spacegroupSubfolders = False):
    fileFormat = fileFormat.lower();

    # Check fileFormat.

    if fileFormat not in SupportedExportFileFormats:
        raise Exception("Error: fileFormat '{0}' is not supported.".format(fileFormat));

    # If workingDirectory does not exist, create it.

    if not os.path.isdir(workingDirectory):
        os.makedirs(workingDirectory);

    # Loop over items in the result set.

    for i, spacegroupGroups in enumerate(resultSet):
        # Sort the spacegropup keys and reorder them to descending symmetry order.

        keys = sorted(spacegroupGroups.keys())[::-1];

        # Determine a chemical formula.
        # We assume here that the supplied resultSet has come from one of the routines in this module, and thus that all structures in each set of spacegroup groups have the same composition.

        structures, _ = spacegroupGroups[keys[0]];

        chemicalFormula = structures[0].GetChemicalFormula(atomicSymbolLookupTable = atomicSymbolLookupTable);

        # Build a name for the archive.

        # If a prefix has been set, start with that.

        archiveName = "{0}_".format(prefix) if prefix != None else "";

        # If there are multiple result sets, include the substitution number in the name.

        if len(resultSet) > 1:
            archiveName = archiveName + "{0:0>3}_".format(i + 1);

        # Finally, append the chemical formula and the file extension.

        archiveName = archiveName + chemicalFormula;

        # Calculate the common divisor to normalise the degeneracies.

        mergedDegeneracies = [];

        for key in keys:
            _, degeneracies = spacegroupGroups[key];
            mergedDegeneracies = mergedDegeneracies + degeneracies;

        commonDivisor = GetCommonDivisor(mergedDegeneracies);

        # Output and archive the structures.

        with tarfile.open(os.path.join(workingDirectory, "{0}.tar.gz".format(archiveName)), 'w:gz') as archiveFile:
            for key in keys:
                spacegroupNumber, spacegroupSymbol = key;

                # If spacegroupSubfolders is set, the structures will be divided into spacegroup subfolders.

                subfolderName = None;

                if spacegroupSubfolders:
                    subfolderName = "{0}-{1}".format(
                        spacegroupNumber, spacegroupSymbol.replace('/', '_')
                        );

                structures, degeneracies = spacegroupGroups[key];

                # Write out each structure and add to the archive.

                for i, (structure, degeneracy) in enumerate(zip(structures, degeneracies)):
                    # Give each structure a title line that includes the chemical formula, spacegroup and normalised degeneracy.

                    structure.SetName(
                        "{0} : SG = {1} ({2}), rel. weight = {3}".format(chemicalFormula, spacegroupNumber, spacegroupSymbol, degeneracy // commonDivisor)
                        );

                    # Generate a file name from the chemical formula, spacegroup number and structure number.

                    fileName = None;

                    if fileFormat == 'vasp':
                        fileName = "{0}_SG-{1}_{2:0>4}.vasp".format(chemicalFormula, spacegroupNumber, i + 1);
                        IO.WritePOSCARFile(structure, fileName);
                    elif fileFormat == 'aims':
                        fileName = "{0}_SG-{1}_{2:0>4}.geometry.in".format(chemicalFormula, spacegroupNumber, i + 1);
                        IO.WriteAIMSGeometryFile(structure, fileName);

                    archiveFile.add(
                        fileName, arcname = "{0}/{1}/{2}".format(archiveName, subfolderName, fileName) if subfolderName != None else "{0}/{1}".format(archiveName, fileName)
                        );

                    # Delete the temporary file once added to the archive.

                    os.remove(fileName);


# ------------------
# Printing Functions
# ------------------

def PrintSpacegroupGroupSummary(spacegroupGroups):
    # Sort the dictionary keys.
    # The first element of the key tuples is the spacegroup number, so sorting will put the keys ascending symmetry order.
    # It's more intuitive to print the table rows in order of descending symmetry, so we reverse the list.

    keys = sorted(spacegroupGroups.keys())[::-1];

    # Obtain the number of structures and the sum of the degeneracies in each group.

    structureCounts, degeneracySums = [], [];

    for key in keys:
        structures, degeneracies = spacegroupGroups[key];

        structureCounts.append(len(structures));
        degeneracySums.append(sum(degeneracies));

    # Work out the maximum integer value to be printed, and hence the required length of the formatted text field.

    maxValue = max(
        max(structureCounts), max(degeneracySums)
        );

    fieldLength = max(
        len("{0:,}".format(maxValue)), 16
        );

    # Print a summary table.

    headerRowFormatCode = "{{0: ^16}} | {{1: ^{0}}} | {{2: ^{0}}}".format(fieldLength);

    headerRow = headerRowFormatCode.format("Spacegroup", "# Structures", "# Unique");

    print(headerRow);
    print('-' * len(headerRow));

    dataRowFormatCode = "{{0: <3}} {{1: <12}} | {{2: >{0},}} | {{3: >{0},}}".format(fieldLength);

    for key, structureCount, degeneracySum in zip(keys, structureCounts, degeneracySums):
        spacegroupNumber, spacegroupSymbol = key;
        print(dataRowFormatCode.format(spacegroupNumber, spacegroupSymbol, degeneracySum,  structureCount));

    print("");


# ----------------------
# Misc Utility Functions
# ----------------------

def GetCommonDivisor(integers):
    commonDivisor = 1;

    hasDivisor = True;

    while hasDivisor:
        # Iteratively divide through by the smallest value until the common divisor found.

        divisor = min(integers);

        if divisor == 1:
            break;

        for value in integers:
            if value % divisor != 0:
                hasDivisor = False;
                break;

        if hasDivisor:
            commonDivisor = commonDivisor * divisor;
            integers = [integer // divisor for integer in integers];

    return commonDivisor;
