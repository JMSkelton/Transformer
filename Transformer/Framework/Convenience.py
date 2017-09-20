# Transformer/Framework/Convenience.py by J. M. Skelton


# -------
# Imports
# -------

import multiprocessing;
import warnings;

from Transformer import Structure;

from Transformer.Framework import Core;

from Transformer.Utilities import StructureTools;


# ----------------
# Antisite Defects
# ----------------

def AntisiteDefects(
        structure, atom1, atom2, numDefects = None,
        tolerance = None,
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

    # Collect structures with antisite defects by performing a sequence of substitutions where we swap atom1 -> 1 and atom2 -> 2.

    _, antisiteDefects = Core.AtomicSubstitutions(
        structure, [(atom1, placeholder1), (atom2, placeholder2)] * numDefects,
        storeIntermediate = [i for i in range(0, 2 * numDefects + 1, 2)], tolerance = tolerance,
        printProgressUpdate = printProgressUpdate,
        **kwargs
        );

    # Modify the substituted structures in the result set to swap the placeholders for 1 -> atom2 and 2 -> atom1.

    for i, spacegroupGroups in enumerate(antisiteDefects[1:]):
        if printProgressUpdate:
            print("AntisiteDefects(): Post processing defect set {0}.".format(i + 1));
            print("");

        # Merge the structures and degeneracies in the spacegroup groups into a "flat" lists.

        structuresFlat, degeneraciesFlat = [], [];

        for structures, degeneracies in spacegroupGroups.values():
            structuresFlat = structuresFlat + structures;
            degeneraciesFlat = degeneraciesFlat + degeneracies;

        # Replace the placeholders.

        for structure in structuresFlat:
            structure.SwapAtoms(placeholder1, atom2);
            structure.SwapAtoms(placeholder2, atom1);

        numStructures = len(structuresFlat);

        # Regroup the structures.

        spacegroupGroups = StructureTools.GroupStructuresBySpacegroup(
            structuresFlat, degeneraciesFlat, tolerance = tolerance
            );

        if printProgressUpdate:
            StructureTools.PrintSpacegroupGroupSummary(spacegroupGroups);

        # Update the results.

        antisiteDefects[i + 1] = spacegroupGroups;

    # Return the result.

    return antisiteDefects;


# ---------------
# Solid Solutions
# ---------------

def SolidSolution(structure, atom1, atom2, useShortcut = True, printProgressUpdate = True, **kwargs):
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

    _, solidSolutions = Core.AtomicSubstitutions(
        structure, substitutions,
        printProgressUpdate = printProgressUpdate,
        **kwargs
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

        for i in range(len(substitutions) + 1, substitutionCount + 1):
            swapIndex = substitutionCount - i;

            # If requested, print a progress update.

            if printProgressUpdate:
                print("SolidSolution(): Inverting substitution set {0} -> {1}".format(swapIndex, i));

            resultSet = solidSolutions[swapIndex];

            newResultSet = { };

            for key, (structures, degeneracies) in resultSet.items():
                newStructures = [];

                for structure in structures:
                    # Clone the structure.

                    structure = structure.Clone();

                    # Get a placeholder to temporarily swap atoms.

                    placeholder = structure.GetAtomTypeNumberPlaceholder();

                    if swapIndex != 0:
                        # Temporarily swap atomReplace with placeholder.

                        structure.SwapAtoms(atomReplace, placeholder);

                    # Swap atomFind with atomReplace.

                    structure.SwapAtoms(atomFind, atomReplace);

                    if swapIndex != 0:
                        # Finally, replace placeholder with atomFind.

                        structure.SwapAtoms(placeholder, atomFind);

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

            # Print a summary of the new results ("repurposing" the _PrintResultSummary() function in the Framework.Core module).

            Core._PrintResultSummary(
                [None] + [(atom1, atom2)] * substitutionCount, solidSolutions, permutationCounts, [i for i in range(0, len(solidSolutions))]
                );

    # Return the result.

    return solidSolutions;
