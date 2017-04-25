# Example_CZTS-AntisiteDefects.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;
from Transformer.ConvenienceFunctions import AntisiteDefects, ExportAtomicSubstitutionResultSet;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("Cu2ZnSnS4.vasp");

# Generate a 2x2x1 supercell (64 atoms).

supercell = structure.GetSupercell((2, 2, 1));

# Generate inequivalent structures with 1 and 2 Cu <-> Zn antisite defects using the AntisiteDefects convenience function.
# At present, this takes 9-10 hrs (!) on a Core i5 iMac (Late 2014).
# This can be _massively_ cut down by setting numDefects = 2 or numDefects = 3.

antisiteDefects = AntisiteDefects(supercell, 'Cu', 'Zn', numDefects = 4);

# Output the results.

ExportAtomicSubstitutionResultSet(antisiteDefects, prefix = "CZTS-AntisiteDefects", workingDirectory = r"Example_CZTS-AntisiteDefects");
