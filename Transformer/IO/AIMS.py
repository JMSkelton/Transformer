# Transformer/IO/AIMS.py by J. M. Skelton


# -------
# Imports
# -------

from Transformer.Structure import Structure;

from Transformer.Utilities import StructureTools;


# ---------
# Functions
# ---------

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

    atomPositions = StructureTools.CartesianToFractionalCoordinates(latticeVectors, atomPositions);

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
