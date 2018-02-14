# Transformer/IO/_VASP.py


# -------
# Imports
# -------

import warnings;

from Transformer import Constants;
from Transformer.IO import _Common;
from Transformer.Utilities import StructureTools;

from Transformer.Structure import Structure;


# ---------
# Functions
# ---------

def ReadPOSCARFile(inputReader, atomTypeNumberLookupTable = None):
    # Variables to collect.

    systemName = None;
    scaleFactor = None;
    latticeVectors = None;
    atomTypes, atomCounts = None, None;
    coordinateType, atomPositions = None, None;

    # Read the system name.

    systemName = next(inputReader).strip();

    # Read the scale factor.

    scaleFactor = float(next(inputReader).strip());

    # Read the lattice vectors.

    latticeVectors = [];

    for i in range(0, 3):
        latticeVectors.append(
            [float(element) for element in next(inputReader).strip().split()][:3]
            );

    # Although we sliced the list returned from split(), this does not guarentee that there were at least three elements.

    for latticeVector in latticeVectors:
        if len(latticeVector) != 3:
            raise Exception("Error: The lattice vector specification in the supplied VASP POSCAR file is invalid.");

    # Read the atom types and/or atom counts.

    atomTypes = [element for element in next(inputReader).strip().split()];

    atomCounts = None;

    if atomTypes[0].isdigit():
        atomCounts = [int(item) for item in atomTypes];
        atomTypes = None;
    else:
        atomCounts = [int(element) for element in next(inputReader).strip().split()];

    # If atom types were given in the file, check the number of atom types listed is consistent with the number of atom counts.

    if atomTypes != None and len(atomTypes) != len(atomCounts):
        raise Exception("Error: The atom-type and atom-count lines in the supplied VASP POSCAR file contain different numbers of entries.");

    # Read the coordinate type.

    coordinateType = None;

    keyword = next(inputReader).strip().lower();

    # Check for and skip the "selective dynamics" keyword.

    if keyword[0] == "s":
        keyword = next(inputReader).strip().lower();

    if keyword[0] == 'd':
        coordinateType = 'd';
    elif keyword[0] == 'c' or keyword[0] == 'k':
        coordinateType = 'c';
    else:
        raise Exception("Error: The coordinate-type line in the supplied VASP POSCAR file contains an unexpected keyword.");

    # Read the atom positions.

    totalAtomCount = 0;

    for atomCount in atomCounts:
        totalAtomCount = totalAtomCount + atomCount;

    atomPositions = [];

    for i in range(0, totalAtomCount):
        elements = next(inputReader).strip().split();

        atomPositions.append(
            [float(element) for element in elements[:3]]
            );

    for atomPosition in atomPositions:
        if len(atomPosition) != 3:
            raise Exception("Error: One or more atomic position specifications in the supplied VASP POSCAR file is invalid.");

    # If a scale factor other than 1 has been set, adjust the lattice vectors.

    if scaleFactor != 1.0:
        for i, vector in enumerate(latticeVectors):
            latticeVectors[i] = [scaleFactor * x for x in vector];

    # Build a list of atom-type numbers.

    atomTypeNumbers = None;

    if atomTypes != None:
        # If atom types were read from the POSCAR file, convert these to atomic numbers.

        atomicSymbols = [];

        for atomType, atomCount in zip(atomTypes, atomCounts):
            atomicSymbols = atomicSymbols + [atomType] * atomCount;

        # Convert the atomic symbols to atom-type numbers.

        atomTypeNumbers = _Common.AtomicSymbolsToAtomTypeNumbers(atomicSymbols, atomTypeNumberLookupTable = atomTypeNumberLookupTable);
    else:
        # If not, issue a warning and assign negative type numbers from -1.

        warnings.warn("Structure objects returned by reading VASP 4-format POSCAR files numbers will be initialised with negative atomic numbers from -1.", UserWarning);

        atomTypeNumbers = [];

        for i, atomCount in enumerate(atomCounts):
            atomTypeNumbers = atomTypeNumbers + [-1 * (i + 1)] * atomCount;

    # If the atom positions are given in Cartesian coordinates, convert them to fractional coordinates.

    if coordinateType == 'c':
        atomPositions = StructureTools.CartesianToFractionalCoordinates(latticeVectors, atomPositions);

    # Return a Structure object.

    return Structure(latticeVectors, atomPositions, atomTypeNumbers, name = systemName);

def WritePOSCARFile(structure, outputWriter, atomicSymbolLookupTable = None):
    # Write the system name; Structure.GetName() returns a sensible default value if a name is not set.

    outputWriter.write("{0}\n".format(structure.GetName()));

    # Write the scale factor.

    outputWriter.write("  {0: >19.16f}\n".format(1.0));

    # Write the lattice vectors.

    for ax, ay, az in structure.GetLatticeVectors():
        outputWriter.write("  {0: >21.16f}  {1: >21.16f}  {2: >21.16f}\n".format(ax, ay, az));

    # Write the atom types and counts.

    atomicSymbols, atomCounts = structure.GetAtomicSymbolsCounts(atomicSymbolLookupTable = atomicSymbolLookupTable);

    for atomicSymbol in atomicSymbols:
        outputWriter.write("  {0: >3}".format(atomicSymbol));

    outputWriter.write("\n");

    for atomCount in atomCounts:
        outputWriter.write("  {0: >3}".format(atomCount));

    outputWriter.write("\n");

    # Write the coordinate type.

    outputWriter.write("Direct\n");

    # Write the atom positions.

    for x, y, z in structure.GetAtomPositions():
        outputWriter.write("  {0: >21.16f}  {1: >21.16f}  {2: >21.16f}\n".format(x, y, z));
