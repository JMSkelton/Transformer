# Transformer/Utilities/StructureTools.py by J. M. Skelton


# -------
# Imports
# -------

import numpy as np;


# ---------
# Functions
# ---------

def CartesianToFractionalCoordinates(latticeVectors, atomPositions):
    # Treat the lattice vectors as a 3x3 matrix and invert to obtain the transformation matrix to fractional coordinates.

    transformationMatrix = np.linalg.inv(latticeVectors);

    # Return the atom positions multiplied by the transformation matrix.

    return [np.dot(position, transformationMatrix) for position in atomPositions];

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
