# Transformer/IO.py by J. M. Skelton


# -------
# Imports
# -------

import warnings;

import numpy as np;

try:
    # spglib version <= 1.8.x

    import spglib as spg;
except ImportError:
    # Nwwer spglib versions.

    from pyspglib import spglib as spg;

from Transformer import Constants;
from Transformer import Utilities;

from Transformer.Structure import Structure;


# -----------------
# VASP Input/Output
# -----------------

def ReadPOSCARFile(filePath):
    # Variables to collect.

    systemName = None;
    scaleFactor = None;
    latticeVectors = None;
    atomTypes, atomCounts = None, None;
    coordinateType, atomPositions = None, None;

    with open(filePath, 'r') as inputReader:
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
                raise Exception("Error: Lattice vector specification in VASP POSCAR file \"{0}\" is invalid.".format(filePath));

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
            raise Exception("Error: The atom-type and atom-count lines in VASP POSCAR file \"{0}\" contain different numbers of entries.".format(filePath));

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
            raise Exception("Error: The coordinate-type line in VASP POSCAR file \"{0}\" contains an unexpected keyword.".format(filePath));

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
                raise Exception("Error: Atomic position specification in VASP POSCAR file \"{0}\" is invalid.".format(filePath));

    # If a scale factor other than 1 has been set, adjust the lattice vectors.

    if scaleFactor != 1.0:
        for i, vector in enumerate(latticeVectors):
            latticeVectors[i] = [scaleFactor * x for x in vector];

    # Build a list of atom-type numbers.

    atomTypeNumbers = [];

    if atomTypes != None:
        # If atom types were read from the POSCAR file, convert these to atomic numbers.

        for atomType, atomCount in zip(atomTypes, atomCounts):
            atomTypeNumbers = atomTypeNumbers + [Constants.SymbolToAtomicNumber(atomType)] * atomCount;
    else:
        # If not, issue a warning and assign sequential type numbers from 1.

        warnings.warn("Structure objects returned by reading VASP 4-format POSCAR files numbers cannot be initialised with atomic numbers and instead use sequential atom-type numbers from 1.", UserWarning);

        for i, atomCount in enumerate(atomCounts):
            atomTypeNumbers = atomTypeNumbers + [i + 1] * atomCount;

    # If the atom positions are given in Cartesian coordinates, convert them to fractional coordinates.

    if coordinateType == 'c':
        atomPositions = Utilities.CartesianToFractionalCoordinates(latticeVectors, atomPositions);

    # Return a Structure object.

    return Structure(
        latticeVectors, atomPositions, atomTypeNumbers, name = systemName
        );

def WritePOSCARFile(structure, filePath, atomicSymbolLookupTable = None):
    with open(filePath, 'w') as outputWriter:
        # Write the system name; Structure.GetName() returns a sensible default value if a name is not set.

        outputWriter.write("{0}\n".format(structure.GetName()));

        # Write the scale factor.

        outputWriter.write("  {0: >19.16f}\n".format(1.0));

        # Write the lattice vectors.

        for ax, ay, az in structure.GetLatticeVectors():
            outputWriter.write("  {0: >21.16f}  {1: >21.16f} {2: >21.16f}\n".format(ax, ay, az));

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
            outputWriter.write("  {0: >21.16f}  {1: >21.16f} {2: >21.16f}\n".format(x, y, z));


# ---------------------
# FHI-AIMS Input/Output
# ---------------------

def ReadAIMSGeometryFile(filePath):
    # Variables to collect.

    latticeVectors = [];
    atomPositions, atomicSymbols = [], [];

    with open(filePath, 'r') as inputReader:
        # Read the file and parse and store lattice_vector and atom keyword entries.

        for line in inputReader:
            # Strip whitespace.

            line = line.strip();

            # Strip comments, if present.

            if '#' in line:
                line = line[:line.find('#')];

            # Test for and parse lattice_vector and atom keywords.

            if len(line) > 14 and line[:14].lower() == "lattice_vector":
                latticeVectors.append(
                    [float(element) for element in line.split()[1:4]]
                    );

            elif len(line) > 4 and line[:4].lower() == "atom":
                elements = line.split();

                atomPositions.append(
                    [float(element) for element in elements[1:4]]
                    );

                atomicSymbols.append(elements[4]);

    # Sanity checks.

    if len(latticeVectors) != 3:
        raise Exception("Error: Incorrect number of lattice_vector keywords in AIMS geometry file \"{0}\".".format(filePath));

    if len(atomPositions) == 0:
        raise Exception("Error: No atom positions/atomic symbols found in AIMS geometry file \"{0}\" - this may be because this function only currently supports the 'atom' keyword.".format(filePath));

    # Convert the atom positions from Cartesian to fractional coordinates.

    atomPositions = Utilities.CartesianToFractionalCoordinates(latticeVectors, atomPositions);

    # Return a Structure object.

    return Structure(
        latticeVectors, atomPositions, atomicSymbols
        );

def WriteAIMSGeometryFile(structure, filePath, atomicSymbolLookupTable = None):
    with open(filePath, 'w') as outputWriter:
        # Write the structure name as a comment.

        outputWriter.write("# {0}\n".format(structure.GetName()));

        # Write the lattice vectors.

        for ax, ay, az in structure.GetLatticeVectors():
            outputWriter.write("lattice_vector {0: >15.8f} {1: >15.8f} {2: >15.8f}\n".format(ax, ay, az));

        # Write the atom positions and types.

        atomicSymbols, atomCounts = structure.GetAtomicSymbolsCounts();

        atomPositions = structure.GetAtomPositionsCartesian();

        atomPositionsPointer = 0;

        for atomicSymbol, atomCount in zip(atomicSymbols, atomCounts):
            for i in range(0, atomCount):
                x, y, z = atomPositions[atomPositionsPointer];
                outputWriter.write("atom  {0: >15.8f} {1: >15.8f} {2: >15.8f} {3}\n".format(x, y, z, atomicSymbol))

                atomPositionsPointer = atomPositionsPointer + 1;
