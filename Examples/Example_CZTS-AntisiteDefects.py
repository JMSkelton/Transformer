# Example_CZTS-AntisiteDefects.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;
from Transformer.ConvenienceFunctions import AntisiteDefects, ExportResultSet;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("Cu2ZnSnS4.vasp");

# Generate a 2x2x1 supercell (64 atoms).

supercell = structure.GetSupercell((2, 2, 1));

# Generate inequivalent structures with 1 and 2 Cu <-> Zn antisite defects using the AntisiteDefects convenience function.

antisiteDefects = AntisiteDefects(supercell, 'Cu', 'Zn', numDefects = 2);

# Output the results.

ExportResultSet(
    antisiteDefects, prefix = "CZTS-AntisiteDefects", workingDirectory = r"Example_CZTS-AntisiteDefects"
    );
