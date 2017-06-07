# Transformer/Structure.py by J. M. Skelton


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


# ---------------
# Structure Class
# ---------------

class Structure:
    # -----------
    # Constructor
    # -----------

    def __init__(self, latticeVectors, atomPositions, atomTypes, name = None):
        # Convert the lattice vectors to a NumPy array, check the shape, and store.

        latticeVectors = np.array(latticeVectors, dtype = np.float64);

        dim1, dim2 = latticeVectors.shape;

        if dim1 != dim2 != 3:
            raise Exception("Error: latticeVectors must be convertible to a 3x3 NumPy array.");

        self._latticeVectors = latticeVectors;

        # Obtain a list of atom-type numbers.

        atomTypeNumbers = [
            AtomTypeToAtomTypeNumber(atomType) for atomType in atomTypes
            ];

        # Check that the length of atomTypeNumbers/atomicSymbols is consistent with that of atomPositions.

        if len(atomTypeNumbers) != len(atomPositions):
            raise Exception("Error: The lengths of atomTypeNumbers/atomicSymbols and atomPositions are not consistent.");

        # Generate a NumPy structured array to store the atom data.
        # The atom positions are adjusted into the range [0, 1], which can be done as a one liner with the non-standrd (?) Python modulus (%) operator.

        atoms = np.array(
            [tuple([typeNumber] + [x % 1.0 for x in position]) for typeNumber, position in zip(atomTypeNumbers, atomPositions)],
            dtype = [('type_number', 'i4'), ('pos_x', 'f8'), ('pos_y', 'f8'), ('pos_z', 'f8')]
            );

        # Sort the data.

        atoms.sort(axis = 0);

        # Store the atoms.

        self._atoms = atoms;

        # Store the name.

        self._name = name;

        # Initialise the fields used to store lazily-initialised properties.

        self._ResetLazyPropertyFields();

    # ---------------
    # Private Methods
    # ---------------

    def _ResetLazyPropertyFields(self):
        # Reset fields used for lazy initialisation of some of the quantities returned by the Get* functions.

        self._pLatticeVectors = None;
        self._pAtomTypeNumbers = None;
        self._pAtomPositions = None;
        self._pAtomPositionsCartesian = None;

        self._pAtomicSymbolsCounts = None;
        self._pChemicalFormula = None;

        self._symmetryAnalysisTolerance = None;

        self._pSpacegroup = None;
        self._pSymmetryOperations = None;
        self._pUniqueAtomIndices = None;

    def _PerformSymmetryAnalysis(self, tolerance = None):
        # If no symmetry tolerance is provided, use the default value.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Only perform the symmetry analysis if it hasn't been done yet (_symmetryAnalysisTolerance == None) or if the tolerance has changed.

        if tolerance != self._symmetryAnalysisTolerance:
            # Call the get_spacegroup() routine from spglib.

            result = spg.get_spacegroup(
                (self.GetLatticeVectors(), self.GetAtomPositions(), self.GetAtomTypeNumbers()),
                symprec = tolerance
                );

            # The result is of the form "<symbol> (<number>)".

            elements = result.split(' ');

            # Exctract the symbol and number, convert the number to an integer, and store both parts as a tuple.

            self._pSpacegroup = (int(elements[-1][1:-1]), ' '.join(elements[:-1]));

            # Call the get_symmetry() routine from spglib.

            result = spg.get_symmetry(
                (self.GetLatticeVectors(), self.GetAtomPositions(), self.GetAtomTypeNumbers()),
                symprec = tolerance
                );

            # If the symmetry search fails, get_symmetry() returns None.

            if result == None:
                raise Exception("Error: spglib get_symmetry() routine returned None.");

            # Store the symmetry operations as a list of tuples.

            self._pSymmetryOperations = [item for item in zip(result['rotations'], result['translations'])];

            # Get the indices of the unique atoms.

            uniqueAtomIndices = np.unique(result['equivalent_atoms']);

            # Store the unique atom indices along with the number of equivalent sites.

            self._pUniqueAtomIndices = (
                uniqueAtomIndices,
                [len(np.where(result['equivalent_atoms'] == index)[0]) for index in uniqueAtomIndices]
                );

            # Store the symmetry tolerance used to call the spglib functions.

            self._symmetryAnalysisTolerance = tolerance;

    # --------------------------
    # Property Get*/Set* Methods
    # --------------------------

    def SetAtom(self, index, atomType, atomPosition = None):
        # Sanity check.

        numAtoms, = self._atoms.shape;

        if index >= numAtoms:
            raise Exception("Error: Index {0} is out of range for number of atoms {1}.".format(index, numAtoms));

        atoms = self._atoms;

        # Convert to atomType number to an atom-type number.

        atomTypeNumber = AtomTypeToAtomTypeNumber(atomType);

        if atomTypeNumber != None:
            # Substitute atom.

            atoms[index]['type_number'] = atomTypeNumber;

            if atomPosition != None:
                atoms[index]['pos_x', 'pos_y', 'pos_z'] = tuple(x % 1.0 for x in atomPosition);

            # Resort atoms.

            atoms.sort();
        else:
            if atomPosition != None:
                warnings.warn("When atomType is set to None, the atom is removed and its position is not updated.", UserWarning);

            # Delete atom.

            atomsNew = np.zeros(numAtoms - 1, dtype = [('type_number', 'i4'), ('pos_x', 'f8'), ('pos_y', 'f8'), ('pos_z', 'f8')]);

            atomsNew[:index] = atoms[:index];
            atomsNew[index:] = atoms[index + 1:];

            # Since atoms should be sorted, removing an entry shouldn't require resorting.

            self._atoms = atomsNew;

        # Invalidate lazily-initialised property fields.

        self._ResetLazyPropertyFields();

    def GetName(self):
        name = self._name;

        if name == None or name == "":
            # If name empty, return a chemical formula if available, or else return the default value.

            if self._pChemicalFormula != None:
                # We don't know a priori whether the atom-type numbers map to elements, or whether these are arbitrarily set.
                # However, if the user has requested a chemical formula without supplying a symbol lookup table, then we can assume this is the case.

                return self._pChemicalFormula;
            else:
                return Structure.DefaultName;
        else:
            return name;

    def SetName(self, name):
        self._name = name;

    def GetLatticeVectors(self):
        # Return the lattice vectors as a list of 1D NumPy arrays.

        if self._pLatticeVectors == None:
            # We use np.copy(), because basic slicing may return a view which refers to the internal lattice vectors.

            self._pLatticeVectors = [np.copy(self._latticeVectors[i, :]) for i in range(0, 3)];

        return self._pLatticeVectors;

    def GetAtomTypeNumbers(self):
        # Return the atom-type numbers as a list of integers.

        if self._pAtomTypeNumbers == None:
            self._pAtomTypeNumbers = [item[0] for item in self._atoms[['type_number']]];

        return self._pAtomTypeNumbers;

    def GetAtomPositions(self):
        # Return the atom positions as a list of 1D NumPy arrays.

        if self._pAtomPositions == None:
            self._pAtomPositions = [
                np.array([item['pos_x'], item['pos_y'], item['pos_z']], dtype = np.float64)
                    for item in self._atoms[['pos_x', 'pos_y', 'pos_z']]
                ];

        return self._pAtomPositions;

    def GetAtomPositionsCartesian(self):
        # Convert the atom positions to Cartesian coordinates and return as a list of NumPy arrays.

        if self._pAtomPositionsCartesian == None:
            # Load the lattice vectors.

            latticeVectors = self._latticeVectors;

            v1, v2, v3 = latticeVectors[0, :], latticeVectors[1, :], latticeVectors[2, :];

            # Load the atom positions and convert to Cartesian coordinates.

            atomPositions = [];

            for item in self._atoms[['pos_x', 'pos_y', 'pos_z']]:
                atomPositions.append(
                    item['pos_x'] * v1 + item['pos_y'] * v2 + item['pos_z'] * v3
                    );

            self._pAtomPositionsCartesian = atomPositions;

        return self._pAtomPositionsCartesian;

    def GetSymmetryAnalysisTolerance(self):
        if self._symmetryAnalysisTolerance == None:
            # If a symmetry tolerance has not been set, return the default value.
            # This is what would be used if the user called one of the Get* functions to obtain symmetry properties without specifying a tolerance.

            return Structure._DefaultSymmetryAnalysisTolerance;
        else:
            return self._symmetryAnalysisTolerance;

    def GetSpacegroup(self, tolerance = None):
        # Perform the symmetry analysis if required.

        self._PerformSymmetryAnalysis(tolerance = tolerance);

        return self._pSpacegroup;

    def GetSymmetryOperations(self, tolerance = None):
        # Perform the symmetry analysis if required.

        self._PerformSymmetryAnalysis(tolerance = tolerance);

        return self._pSymmetryOperations;

    def GetUniqueAtomIndices(self, tolerance = None):
        # Perform the symmetry analysis if required.

        self._PerformSymmetryAnalysis(tolerance = tolerance);

        return self._pUniqueAtomIndices;

    # ------------
    # Core Methods
    # ------------

    def Clone(self):
        # Given the conversions in the constructor, this should produce a deep copy -- at least of the lattice vectors, atom positions and atom-type numbers.

        return Structure(
            self.GetLatticeVectors(), self.GetAtomPositions(), self.GetAtomTypeNumbers(), name = self._name
            );

    def CompareLatticeVectors(self, structure, tolerance = None):
        # If no equivalence tolerance is supplied, use the default.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        return ((structure._latticeVectors - self._latticeVectors) <= tolerance).all();

    def CompareAtomTypeNumbers(self, structure):
        # First check that the supplied structure has the same number of atoms.

        numAtoms1, = structure._atoms.shape;
        numAtoms2, = self._atoms.shape;

        if numAtoms1 != numAtoms2:
            return False;

        # Get a view to the 'type_number' fields in both structures -- this allows them to be compared using NumPy routines.

        types1 = self._atoms[:][['type_number']].view(np.int32);
        types2 = structure._atoms[:][['type_number']].view(np.int32);

        return (types2 == types1).all();

    def CompareAtomPositions(self, structure, tolerance = None):
        # If no equivalence tolerance is supplied, use the default.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Check the supplied structure has the same number of atoms.

        numAtoms1, = structure._atoms.shape;
        numAtoms2, = self._atoms.shape;

        if numAtoms1 != numAtoms2:
            return False;

        # Generate views to the 'pos_*' fields in both structures, and compare using NumPy.

        positions1 = self._atoms[:][['pos_x', 'pos_y', 'pos_z']].view(dtype = np.float64);
        positions2 = structure._atoms[:][['pos_x', 'pos_y', 'pos_z']].view(dtype = np.float64);

        return ((positions2 - positions1) <= tolerance).all();

    def CompareCell(self, structure, tolerance = None):
        # If no equivalence tolerance is supplied, use the default.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Fastest check: the lattice vectors are the same.

        if not self.CompareLatticeVectors(structure, tolerance = tolerance):
            return False;

        # Next fastest check: compare the atom counts and atom-type numbers.

        if not self.CompareAtomTypeNumbers(structure):
            return False;

        # Finally, compare the atom positions.

        if not self.CompareAtomPositions(structure, tolerance = tolerance):
            return False;

        # All tests passed -> the cells are equivalent.

        return True;

    def Compare(self, structure, tolerance = None):
        # Compare names.

        if self._name != structure._name:
            return False;

        # Compare cells.

        if not self.CompareCell(structure, tolerance):
            return False;

        # Tests passed -> Structure objects hold equivalent data.

        return True;

    def GetAtomicSymbolsCounts(self, atomicSymbolLookupTable = None):
        # If atomicSymbolLookupTable is not provided and the internal _pAtomSymbolsCounts has been initialised, we can simply return a copy of that.

        if atomicSymbolLookupTable == None and self._pAtomicSymbolsCounts != None:
            atomicSymbols, atomCounts = self._pAtomicSymbolsCounts;

            return (
                [symbol for symbol in atomicSymbols],
                [atomCount for atomCount in atomCounts]
                );

        atomicSymbolsList = [];

        # If an atomic-symbol lookup table is provided, get the symbols from that; if not get them from the periodic table in the Constants module.

        if atomicSymbolLookupTable != None:
            atomicSymbolsList = [atomicSymbolLookupTable[typeNumber] for typeNumber in self.GetAtomTypeNumbers()];
        else:
            atomicSymbolsList = [Constants.AtomicNumberToSymbol(typeNumber) for typeNumber in self.GetAtomTypeNumbers()];

        # Use the list of types to build a list of atomicSymbols and counts to add to the POSCAR header.

        atomicSymbols, atomCounts = [atomicSymbolsList[0]], [1];

        # The atom-type numbers returned by GetAtomTypeNumbers() should be sorted.

        for i, atomicSymbol in enumerate(atomicSymbolsList[1:]):
            if atomicSymbol == atomicSymbols[-1]:
                atomCounts[-1] = atomCounts[-1] + 1;
            else:
                atomicSymbols.append(atomicSymbol);
                atomCounts.append(1);

        # To avoid having to store a reference to the symbol lookup table, a copy of the atom types/counts are only stored internally if atomicSymbolLookupTable is not set.

        if atomicSymbolLookupTable == None:
            self._pAtomicSymbolsCounts = (
                [symbol for symbol in atomicSymbols],
                [atomCount for atomCount in atomCounts]
                );

        return (atomicSymbols, atomCounts);

    def GetChemicalFormula(self, atomicSymbolLookupTable = None):
        # If atomicSymbolLookupTable is not provided and the internal _pChemicalFormula has been initialised, return a copy of that.

        if atomicSymbolLookupTable == None and self._pChemicalFormula != None:
            return self._pChemicalFormula;

        atomicSymbols, atomCounts = self.GetAtomicSymbolsCounts(atomicSymbolLookupTable = atomicSymbolLookupTable);

        # Obligatory cryptic one-liner (!).

        chemicalFormula = "".join("{0}{1}".format(*item) for item in zip(atomicSymbols, atomCounts));

        # Only store a copy of the chemical formula if atomicSymbolLookupTable is not set.

        if atomicSymbolLookupTable == None:
            self._pChemicalFormula = chemicalFormula;

        return self._pChemicalFormula;

    # ----------------------
    # Transformation Methods
    # ----------------------

    def GetSymmetryTransform(self, rotation, translation):
        # Get and transform the atom positions.

        newAtomPositions = [
            np.dot(rotation, position) + translation for position in self.GetAtomPositions()
            ];

        # Return a new Structure object with the new atom positions.

        return Structure(
            self.GetLatticeVectors(), newAtomPositions, self.GetAtomTypeNumbers()
            );

    def GetSupercell(self, supercellDim):
        dimA, dimB, dimC = supercellDim;

        # Sanity check.

        if dimA < 1 or dimB < 1 or dimC < 1:
            raise Exception("Error: Supercell dimensions must be >= 1.");

        if dimA % 1 != 0 or dimB % 1 != 0 or dimC % 1 != 0:
            raise Exception("Error: Supercell dimensions must be integers.");

        # Check whether we actually need to do anything.

        if dimA > 1 or dimB > 1 or dimC > 1:
            atomPositions = self.GetAtomPositions();
            atomTypeNumbers = self.GetAtomTypeNumbers();

            # Rescale the atom positions.

            divisor = np.array((dimA, dimB, dimC), dtype = np.float64);

            atomPositions = [position / divisor for position in atomPositions];

            # Generate the new atom positions and list of atomic numbers for the supercell.

            newAtomPositions, newAtomTypeNumbers = [], [];

            for i in range(0, dimA):
                for j in range(0, dimB):
                    for k in range(0, dimC):
                        translationVector = np.array(
                            (float(i) / dimA, float(j) / dimB, float(k) / dimC), dtype = np.float64
                            );

                        newAtomPositions = newAtomPositions + [position + translationVector for position in atomPositions];
                        newAtomTypeNumbers = newAtomTypeNumbers + atomTypeNumbers;

            # Generate the new lattice vectors.

            latticeVectors = self._latticeVectors;

            newLatticeVectors = [
                latticeVectors[0, :] * dimA,
                latticeVectors[1, :] * dimB,
                latticeVectors[2, :] * dimC
                ];

            # Generate a new name..

            newName = None;

            if self._name != None:
                # If the current structure has a name, generate a new one by appending "(<dimA>x<dimB>x<dimC> SC)" to it.

                newName = "{0} ({1}x{2}x{3} SC)".format(self._name, dimA, dimB, dimC);

            # Return a new structure.

            return Structure(
                newLatticeVectors, newAtomPositions, newAtomTypeNumbers, name = newName
                );
        else:
            # If dimA = dimB = dimC = 1, simply return a clone.

            return self.Clone();

    def GetAtomSwap(self, atomType1, atomType2):
        # Convert atomType1 and atomType2 to atom-type numbers.

        atomTypeNumber1 = AtomTypeToAtomTypeNumber(atomType1);
        atomTypeNumber2 = AtomTypeToAtomTypeNumber(atomType2);

        # Get the atom-type numbers and collect the indices of atoms to swap.

        atomTypeNumbers = self.GetAtomTypeNumbers();

        swapIndices = [
            i for i, atomTypeNumber in enumerate(atomTypeNumbers)
                if atomTypeNumber == atomTypeNumber1
            ];

        if len(swapIndices) == 0:
            raise Exception("Error: The atom to swap specified by atomType1 was not found in the current structure.");

        if atomTypeNumber2 != None:
            # If atomTyoeNumber2 was supplied via atomTypeNumber2/atomicSymbol2, adjust the atom-type numbers and return a new structure.

            if atomTypeNumber2 == atomTypeNumber1:
                # If the two atom-type numbers are the same, we don't need to do anything -> return a clone.

                return self.Clone();
            else:
                # Make a deep copy of the atom-type numbers -- changing the atom-type numbers in the list returned by GetAtomPositions() will modify the internal cached list.

                newAtomTypeNumbers = [typeNumber for typeNumber in atomTypeNumbers];

                for index in swapIndices:
                    newAtomTypeNumbers[index] = atomTypeNumber2;

                # Return a new structure with the updated atom-type numbers.

                return Structure(
                    self.GetLatticeVectors(), self.GetAtomPositions(), newAtomTypeNumbers
                    );
        else:
            # If atomTypeNumber2 is None, return a new structure with the atoms removed.

            atomPositions = self.GetAtomPositions();

            return Structure(
                self.GetLatticeVectors(),
                [position for i, position in enumerate(atomPositions) if i not in swapIndices],
                [typeNumber for i, typeNumber in enumerate(atomTypeNumbers) if i not in swapIndices]
                );

    # -------------
    # Static Fields
    # -------------

    DefaultName = "Unknown Structure";
    DefaultSymmetryEquivalenceTolerance = 1.0e-5;


# --------------
# Static Methods
# --------------

# This should ideally be attached to the Structure class, but Python 2.x doesn't allow classes to have static methods.

def AtomTypeToAtomTypeNumber(atomType):
    if atomType != None:
        # If atomType is not None, try and convert it to an atom-type number.

        try:
            # Assume atomType is (convertible to) an integer.

            return int(atomType);
        except ValueError:
            # If that fails, assume atomType is an atomic symbol and lookup the corresponding atomic number.

            return Constants.SymbolToAtomicNumber(str(atomType));
    else:
        return None;
