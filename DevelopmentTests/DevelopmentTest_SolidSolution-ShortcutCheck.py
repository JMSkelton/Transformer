# DevelopmentTest_SolidSolution-ShortcutCheck.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;

from Transformer.ConvenienceFunctions.BatchIO import ExportResultSet, ImportResultSet;
from Transformer.ConvenienceFunctions.Substitutions import SolidSolution;

from Transformer.DevelopmentTools import MapResultSetStructures;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("../Examples/SnS-Pnma.vasp");

# Generate a 2x1x2 supercell (32 atoms).

supercell = structure.GetSupercell((2, 1, 2));

# Evaluate the full set of S/Se solid soliutions with useShortcut = False and export the results.

solidSolutions = SolidSolution(supercell, 'S', 'Se', useShortcut = False);

ExportResultSet(
    solidSolutions, prefix = "SnS-Se-SolidSolution", workingDirectory = r"DevelopmentTest_SolidSolution-ShortcutCheck"
    );

# Evaluate the solid soliutions again with useShortcut = True and export.

solidSolutions = SolidSolution(supercell, 'S', 'Se', useShortcut = True);

ExportResultSet(
    solidSolutions, prefix = "SnS-Se-SolidSolution-Shortcut", workingDirectory = r"DevelopmentTest_SolidSolution-ShortcutCheck"
    );

# Read in both sets of results.

resultSets1 = ImportResultSet(prefix = "SnS-Se-SolidSolution", directory = r"DevelopmentTest_SolidSolution-ShortcutCheck");
resultSets2 = ImportResultSet(prefix = "SnS-Se-SolidSolution-Shortcut", directory = r"DevelopmentTest_SolidSolution-ShortcutCheck");

# Run an exhaustive setwise comparison.

for i, (resultSet1, resultSet2) in enumerate(zip(resultSets1, resultSets2)):
    print("Comparing result set {0}...:".format(i + 1));

    structures1, degeneracies1 = resultSet1;
    structures2, degeneracies2 = resultSet2;

    # Print the structure counts.

    print("  -> INFO: Set 1 contains {0} structure(s)".format(len(structures1)));
    print("  -> INFO: Set 2 contains {0} structure(s)".format(len(structures2)));

    warnings = False;

    # Check the two sets contain the same numbers of structures.

    if len(structures1) != len(structures2):
        print("  -> WARNING: Set 1 and Set 1 contain different numbers of structures.");

        warnings = True;

    # Perform the mapping.

    structureMappings = MapResultSetStructures(supercell, structures1, structures2);

    # Check each structure in the first set maps to one in the second set.

    for j, structureMapping in enumerate(structureMappings):
        numMappings = len(structureMapping);

        if numMappings != 1:
            print("  -> WARNING: Set 1/Structure {0} maps to {1} structure(s)".format(j + 1, numMappings));

            warnings = True;

    # Check each structure in the second set is included in the mapped by something from the first.

    indices2 = set();

    for structureMapping in structureMappings:
        indices2.update(structureMapping);

    for j in range(0, len(structures2)):
        if j not in indices2:
            print("  -> WARNING: Set 2/Structure {0} does not map to anything".format(j + 1));

            warnings = True;

    # Where we have 1:1 mapping, compare the degeneracies.

    for j, structureMapping in enumerate(structureMappings):
        if len(structureMapping) == 1:
            mapIndex = structureMapping[0];

            degeneracy1, degeneracy2 = degeneracies1[j], degeneracies2[mapIndex];

            if degeneracy1 != degeneracy2:
                print("  -> WARNING: Set1/Structure {0} and Set2/Structure {1} degeneracies do not match".format(j + 1, mapIndex + 1));

                warnings = True;

    # If there were no warnings, print an "OK" message.

    if not warnings:
        print("  -> OK: 1:1 mapping between structure(s) and degeneracy(s)");

    print("");
