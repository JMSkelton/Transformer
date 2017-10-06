# Example_SnS-Se-SolidSolution.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO.VASP import ReadPOSCARFile;

from Transformer.Framework.BatchIO import ExportResultSet;
from Transformer.Framework.Convenience import SolidSolution;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("SnS-Pnma.vasp");

# Generate a 2x1x2 supercell (32 atoms).

supercell = structure.GetSupercell((2, 1, 2));

# Evaluate the full set of S/Se solid soliutions.
# This takes < 1 min on a Core i5 iMac (Late 2014) with Cython installed to allow the optimised routines to be used.

solidSolutions = SolidSolution(supercell, 'S', 'Se');

# Output the results.

ExportResultSet(
    solidSolutions, prefix = "SnS-Se-SolidSolution", workingDirectory = r"Example_SnS-Se-SolidSolution"
    );
