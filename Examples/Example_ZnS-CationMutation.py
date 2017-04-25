# Example_ZnS-CationMutation.py by J. M. Skelton

# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;
from Transformer.ConvenienceFunctions import AtomicSubstitutions, PrintSpacegroupGroupSummary, ExportAtomicSubstitutionResultSet;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("ZnS-Conventional.vasp");

# Generate a 1x1x2 supercell (16 atoms).

supercell = structure.GetSupercell((1, 1, 2));

# Define a sequence of substitutions.

substitutions = [];

# First, we replace half of the Zn atoms with Cu, and half with Ga.

substitutions = substitutions + [('Zn', 'Cu')] * 4;
substitutions = substitutions + [('Zn', 'Ga')] * 4;

# Second, replace half the Ga with Zn , and half with Sn.

substitutions = substitutions + [('Ga', 'Zn')] * 2;
substitutions = substitutions + [('Ga', 'Sn')] * 2;

# Use the AtomicSubstitutions convenience function to perform the substitution.
# This should take a very short time to run.
# The storeIntermediate parameter specifies which sets of intermediate structures we want to be returned in the result set; in this case, we ask for the initial structure (0), and the results after the first and second sets of substitutions (4 and 8, respectively).

_, resultSet = AtomicSubstitutions(
    supercell, substitutions, storeIntermediate = [0, 8, 12]
    );

# Print a summary of each set of results.

print("Initial structure:");
print("");

PrintSpacegroupGroupSummary(resultSet[0]);

print("After first substitution:");
print("");

# This summary table should show one unique I-42d structure, which is listed in a COD search for Cu-Ga-S.

PrintSpacegroupGroupSummary(resultSet[1]);

print("After second substitution:");
print("");

# This summary table should show one unique I-42m and one I-4 structure, which correspond to the Stannite and Kesterite structures of CZTS; are both listed in a COD search for Cu-Zn-Sn-S.

PrintSpacegroupGroupSummary(resultSet[2]);

# Export the results.

ExportAtomicSubstitutionResultSet(resultSet, prefix = "ZnS-CationMutation", workingDirectory = r"Example_ZnS-CationMutation");
