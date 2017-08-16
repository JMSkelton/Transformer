# Example_SnS-Se-SolidSolution.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;

from Transformer.ConvenienceFunctions.BatchIO import ExportResultSet;
from Transformer.ConvenienceFunctions.Substitutions import SolidSolution;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("SnS-Pnma.vasp");

# Generate a 2x1x2 supercell (32 atoms).

supercell = structure.GetSupercell((2, 1, 2));

# Evaluate the full set of S/Se solid soliutions.
# This takes 1-2 mins on a Core i5 iMac (Late 2014).

solidSolutions = SolidSolution(supercell, 'S', 'Se');

# Output the results.

ExportResultSet(
    solidSolutions, prefix = "SnS-Se-SolidSolution", workingDirectory = r"Example_SnS-Se-SolidSolution"
    );
