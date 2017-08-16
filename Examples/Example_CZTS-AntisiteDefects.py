# Example_CZTS-AntisiteDefects.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;

from Transformer.ConvenienceFunctions.BatchIO import ExportResultSet;
from Transformer.ConvenienceFunctions.Substitutions import AntisiteDefects;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("Cu2ZnSnS4.vasp");

# Generate a 2x2x1 supercell (64 atoms).

supercell = structure.GetSupercell((2, 2, 1));

# Generate inequivalent structures with up to four Cu <-> Zn antisite defects using the AntisiteDefects convenience function.
# Note: This example is quite heavy, and may take a long time to run without Cython installed to allow the Cython-optimised routines to be used.

antisiteDefects = AntisiteDefects(
    supercell, 'Cu', 'Zn', numDefects = 4
    );

# Output the results.

ExportResultSet(
    antisiteDefects, prefix = "CZTS-AntisiteDefects", workingDirectory = r"Example_CZTS-AntisiteDefects"
    );
