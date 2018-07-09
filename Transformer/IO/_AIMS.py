# Transformer/IO/_AIMS.py


# -------
# Imports
# -------

from Transformer.Utilities import StructureTools;

from Transformer.Structure import Structure;



# ---------
# Functions
# ---------

def ReadGeometryInFile(inputReader, atomicSymbolLookupTable = None):
    # Variables to collect.

    latticeVectors = [];
    atomPositions, atomicSymbols = [], [];

    # Read and parse the file.

    atomPositionsTemp = [];

    for line in inputReader:
        # Strip whitespace.

        line = line.strip();

        # Strip comments, if present.

        if '#' in line:
            line = line[:line.find('#')];

        # Test for and parse lattice_vector and atom[_frac] keywords.

        elements = line.split();

        if len(elements) > 0:
            keyword = elements[0].lower();

            if keyword == "lattice_vector":
                latticeVectors.append(
                    [float(element) for element in elements[1:4]]
                    );

            elif keyword == "atom" or keyword == "atom_frac":
                atomPosition = [float(element) for element in elements[1:4]];

                # Temporarily record the atom positions along with a tag indicating whether they are in fractional coordinates.
                # This is because we may not have the lattice vectors to convert Cartesian positions until we have finished parsing.

                atomPositionsTemp.append(
                    (atomPosition, keyword == "atom_frac")
                    );

                atomicSymbols.append(elements[4]);

    # Sanity checks.

    if len(latticeVectors) != 3:
        raise Exception("Error: An incorrect number of lattice_vector keywords were found in the supplied AIMS geometry file");

    if len(atomPositionsTemp) == 0:
        raise Exception("Error: No atom positions/atomic symbols found in the supplied AIMS geometry file.");

    # Convert any atom positions specified in Cartesian coordinates to fractional coordinates.

    for atomPosition, isFractional in atomPositionsTemp:
        if not isFractional:
            atomPosition = StructureTools.CartesianToFractionalCoordinates(latticeVectors, [atomPosition])[0];

        atomPositions.append(atomPosition);

    # Return a Structure object.

    return Structure(latticeVectors, atomPositions, atomicSymbols, atomicSymbolLookupTable = atomicSymbolLookupTable);

def WriteGeometryInFile(structure, outputWriter, atomicSymbolLookupTable = None):
    # Write the structure name as a comment.

    outputWriter.write("# {0}\n".format(structure.GetName()));

    # Write the lattice vectors.

    for ax, ay, az in structure.GetLatticeVectors():
        outputWriter.write("lattice_vector {0: >15.8f} {1: >15.8f} {2: >15.8f}\n".format(ax, ay, az));

    # Write the atom positions and types.

    atomicSymbols, atomCounts = structure.GetAtomicSymbolsCounts(atomicSymbolLookupTable = atomicSymbolLookupTable);

    atomPositions = structure.GetAtomPositionsCartesian();

    atomPositionsPointer = 0;

    for atomicSymbol, atomCount in zip(atomicSymbols, atomCounts):
        for i in range(0, atomCount):
            x, y, z = atomPositions[atomPositionsPointer];
            outputWriter.write("atom  {0: >15.8f} {1: >15.8f} {2: >15.8f} {3}\n".format(x, y, z, atomicSymbol))

            atomPositionsPointer = atomPositionsPointer + 1;
