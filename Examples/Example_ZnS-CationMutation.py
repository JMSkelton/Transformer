# Example_ZnS-CationMutation.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO.StructureIO import ReadStructure;

from Transformer.Framework.BatchIO import ExportResultSet;
from Transformer.Framework.Core import AtomicSubstitutions;

from Transformer.Utilities.StructureTools import PrintStructureSetSummary;


# ----
# Main
# ----

# Read the input structure.

structure = ReadStructure("ZnS-Conventional.vasp");

# Generate a 1x1x2 supercell (16 atoms).

supercell = structure.GetSupercell(1, 1, 2);

# Define a list of substitution sequences.

substitutions = [];

# First, we replace half of the Zn atoms with Cu, and half with Ga.

substitutions.append(
    [('Zn', 'Cu')] * 4 + [('Zn', 'Ga')] * 4
    );

# Second, replace half the Ga with Zn , and half with Sn.

substitutions.append(
    [('Ga', 'Zn')] * 2 + [('Ga', 'Sn')] * 2
    );

# Use the AtomicSubstitutions convenience function to perform the substitution.
# This should take a very short time to run.

resultSet = AtomicSubstitutions(
    supercell, substitutions
    );

# Print a summary of each set of results.

print("Initial structure:");
print("");

PrintStructureSetSummary(resultSet[0]);

print("After first substitution:");
print("");

# This summary table should show one unique I-42d structure, which is listed in a COD search for Cu-Ga-S.

PrintStructureSetSummary(resultSet[1]);

print("After second substitution:");
print("");

# This summary table should show one unique I-42m and one I-4 structure, which correspond to the Stannite and Kesterite structures of CZTS; are both listed in a COD search for Cu-Zn-Sn-S.

PrintStructureSetSummary(resultSet[2]);

# Export the results.

ExportResultSet(
    resultSet, prefix = "ZnS-CationMutation", archiveDirectory = r"Example_ZnS-CationMutation"
    );
