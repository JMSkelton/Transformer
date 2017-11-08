# Transformer/Framework/Convenience.py


# -------
# Imports
# -------

import multiprocessing;
import warnings;

from Transformer import Structure;
from Transformer import StructureSet;

from Transformer.Framework import Core;

from Transformer.Utilities import StructureTools;


# -------------------------
# Antisite Defects Function
# -------------------------

def AntisiteDefects(
    structure, atom1, atom2, numDefects = None,
    printProgressUpdate = True,
    **kwargs
    ):

    # Set up a generator function and yield the structure sets.

    generator = _AntisiteDefectsIter(
        structure, atom1, atom2, numDefects = numDefects,
        printProgressUpdate = printProgressUpdate,
        **kwargs
        );

    # Run the generator and capture the output.

    substitutions, structureSets, permutationCounts = [], [], [];

    for substitution, structureSet, permutationCount in generator:
        substitutions.append(substitution);

        if structureSet != None:
            structureSets.append(structureSet);

        permutationCounts.append(permutationCount);

    # If required, print a final summary of the results.

    if printProgressUpdate:
        Core._PrintResultSummary(
            substitutions, structureSets, permutationCounts, [i for i in range(0, len(substitutions), 2)]
            );

    return structureSets;

def AntisiteDefectsIter(
    structure, atom1, atom2, numDefects = None,
    printProgressUpdate = True,
    **kwargs
    ):

    # Set up a generator function and yield the structure sets.

    generator = _AntisiteDefectsIter(
        structure, atom1, atom2, numDefects = numDefects,
        printProgressUpdate = printProgressUpdate,
        **kwargs
        );

    for _, structureSet, _ in generator:
        if structureSet != None:
            yield structureSet;

def _AntisiteDefectsIter(
    structure, atom1, atom2, numDefects = None,
    printProgressUpdate = True,
    **kwargs
    ):

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

    placeholder1 = structure.GetAtomTypeNumberPlaceholder();
    placeholder2 = placeholder1 - 1;

    # Build a generator for performing sequences of substitutions where we swap atom1 -> 1 and atom2 -> 2.

    generator = Core._AtomicSubstitutionsIter(
        structure, [[(atom1, placeholder1), (atom2, placeholder2)]] * numDefects,
        printProgressUpdate = printProgressUpdate,
        **kwargs
        );

    # Run the generator.

    counter = 0;

    for substitution, structureSet, permutationCount in generator:
        # For each structure set returned, modify the substituted structures to swap the placeholders for 1 -> atom2 and 2 -> atom1.

        if structureSet != None:
            if printProgressUpdate:
                print("AntisiteDefects: Post processing structure set {0}.".format(counter + 1));
                print("");

            # Merge the structures and degeneracies in the spacegroup groups into a "flat" lists.

            structuresFlat, degeneraciesFlat = structureSet.GetStructureSetFlat();

            # Clone the structures and replace the placeholders.

            newStructures = [];

            for structure in structuresFlat:
                structure = structure.Clone();

                structure.SwapAtoms(
                    [placeholder1, placeholder2], [atom2, atom1]
                    );

                newStructures.append(structure);

            # Put the modified structures into a new structure set.

            newStructureSet = structureSet.CloneNew(
                structures = newStructures, degeneracies = degeneraciesFlat, noInitialMerge = True
                );

            if printProgressUpdate:
                StructureTools.PrintStructureSetSummary(newStructureSet);

            structureSet = newStructureSet;

            counter = counter + 1;

        # Yield the result.

        yield (substitution, structureSet, permutationCount);


# -----------------------
# Solid Solution Function
# -----------------------

def SolidSolution(structure, atom1, atom2, printProgressUpdate = True, **kwargs):
    # Set up generator function.

    generator = _SolidSolutionIter(
        structure, atom1, atom2, printProgressUpdate = printProgressUpdate, **kwargs
        );

    # Run the generator and capture the output.

    substitutions, structureSets, permutationCounts = [], [], [];

    for substitution, structureSet, permutationCount in generator:
        substitutions.append(substitution);
        structureSets.append(structureSet);
        permutationCounts.append(permutationCount);

    # If required, print a final summary of the results.

    if printProgressUpdate:
        Core._PrintResultSummary(
            substitutions, structureSets, permutationCounts, [i for i in range(0, len(structureSets))]
            );

    return structureSets;

def SolidSolutionIter(structure, atom1, atom2, printProgressUpdate = True, **kwargs):
    # Set up a generator function.

    generator = _SolidSolutionIter(
        structure, atom1, atom2, printProgressUpdate = printProgressUpdate, **kwargs
        );

    # Run the generator and yield the structure sets.

    for _, structureSet, _ in generator:
        yield structureSet;

def _SolidSolutionIter(structure, atom1, atom2, printProgressUpdate = True, **kwargs):
    # Convert atom1 and atom2 to atom-type numbers.

    atomFind = Structure.AtomTypeToAtomTypeNumber(atom1);
    atomReplace = Structure.AtomTypeToAtomTypeNumber(atom2);

    atomTypeNumbers = structure.GetAtomTypeNumbers();

    if atomFind not in atomTypeNumbers:
        raise Exception("Error: The atom specified by atomTypeNumber1/atomicSymbol1 was not found in the supplied structure.");

    if atomReplace in atomTypeNumbers:
        raise Exception("Error: The atom specified by atomTypeNumber2/atomicSymbol2 was found in the supplied structure - this is most likely an error.");

    # Count atoms to substitute.

    substitutionCount = 0;

    for atomType in atomTypeNumbers:
        if atomType == atomFind:
            substitutionCount = substitutionCount + 1;

    # We only need to enumerate the structures for up to ~50 % substitution, then we can generate the rest by swapping atoms.

    substitutions = substitutions = [(atom1, atom2)] * (substitutionCount // 2);

    # Perform substitutions.

    solidSolutions, permutationCounts = [], [];

    generator = Core._AtomicSubstitutionsIter(
        structure, substitutions,
        printProgressUpdate = printProgressUpdate,
        **kwargs
        );

    for _, structureSet, permutationCount in generator:
        solidSolutions.append(structureSet);
        permutationCounts.append(permutationCount);

        # Yield the substitution, structure set and permutation count.

        yield ((atom1, atom2), structureSet, permutationCount);

    # Generate the structures for the remaining substitutions.

    for i in range(len(substitutions) + 1, substitutionCount + 1):
        swapIndex = substitutionCount - i;

        if printProgressUpdate:
            print("SolidSolution: Inverting substitution set {0} -> {1}".format(swapIndex, i));

        structureSet = solidSolutions[swapIndex];
        structuresFlat, degeneraciesFlat = structureSet.GetStructureSetFlat();

        newStructures = [];

        for structure, degeneracy in zip(structuresFlat, degeneraciesFlat):
            # Clone the structure, swap atomFind and atomReplace, then add it to the new list of structures.

            structure = structure.Clone();

            structure.SwapAtoms(
                [atomFind, atomReplace], [atomReplace, atomFind]
                );

            # Add the new structure to the list.

            newStructures.append(structure);

        # Build and yield the new structure set.

        newStructureSet = structureSet.CloneNew(
            structures = newStructures, degeneracies = degeneraciesFlat, noInitialMerge = True
            );

        solidSolutions.append(newStructureSet);

        yield ((atom1, atom2), newStructureSet, permutationCounts[swapIndex]);

    print("");
