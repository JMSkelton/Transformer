# Transformer/Structure.py by J. M. Skelton


# -------
# Imports
# -------

import math;
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

        numAtoms = len(atomTypeNumbers);

        # The atom data is stored in a NumPy structured array, which enables NumPy sorting.

        atomData = np.zeros(numAtoms, dtype = Structure._AtomDataType);

        # Creating views to the atom-type numbers and atom positions allows easier manipulation of the data.

        atomTypeNumbersView = atomData.view(dtype = np.int64).reshape((numAtoms, 4))[:, 0];

        atomPositionsView = atomData.view(dtype = np.float64).reshape((numAtoms, 4))[:, 1:];

        # Set the data, clamp the atom positions to the range [0, 1] and sort.

        atomTypeNumbersView[:] = atomTypeNumbers;

        atomPositionsView[:] = atomPositions;
        atomPositionsView[:] %= 1.0;

        atomData.sort();

        # Store the atom data and views.

        self._atomData = atomData;

        self._atomTypeNumbersView = atomTypeNumbersView;
        self._atomPositionsView = atomPositionsView;

        # Store the name.

        self._name = name;

        # Initialise the fields used to store lazily-initialised properties.

        self._ResetLazyPropertyFields();

    # ---------------
    # Private Methods
    # ---------------

    def _ResetLazyPropertyFields(self):
        # Reset fields used for lazy initialisation of some quantities that are non-trivial to compute.

        self._symmetryAnalysisTolerance = None;

        self._pSpacegroup = None;
        self._pSymmetryOperations = None;
        self._pUniqueAtomIndices = None;

        self._pAtomPositionsCartesian = None;

        self._pAtomicSymbolsCounts = None;
        self._pChemicalFormula = None;

    def _PerformSymmetryAnalysis(self, tolerance = None):
        # If no symmetry tolerance is provided, use the default value.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Only perform the symmetry analysis if it hasn't been done yet (_symmetryAnalysisTolerance == None) or if the tolerance has changed.

        if tolerance != self._symmetryAnalysisTolerance:
            # Call the get_spacegroup() routine from spglib.

            result = spg.get_symmetry_dataset(
                (self._latticeVectors, self._atomPositionsView, self._atomTypeNumbersView),
                symprec = tolerance
                );

            # If the symmetry search fails, get_symmetry_dataset() returns None.

            if result == None:
                raise Exception("Error: spglib get_symmetry_dataset() routine returned None.");

            # Stpre the spacegroup number and international symbol..

            self._pSpacegroup = (result['number'], result['international']);

            # spglib returns symmetry operations in the form of a rotation matrix plus a translation.
            # The translation vectors sometimes show numerical noise, which can sometimes cause structures related by symmetry to be detected as inequivalent.
            # We therefore round them the translations based on the specified tolerance and clamp to the range [0, 1].

            roundDecimals = -1 * int(
                math.floor(math.log10(tolerance))
                );

            translations = [
                np.round(translation, roundDecimals) % 1.0
                    for translation in result['translations']
                ];

            # Store the symmetry operations as a list of (rotation, translation) tuples.

            self._pSymmetryOperations = [item for item in zip(result['rotations'], translations)];

            # Get the indices of the unique atoms.

            uniqueAtomIndices = np.unique(result['equivalent_atoms']);

            # Store the unique atom indices along with the number of equivalent sites.

            self._pUniqueAtomIndices = (
                uniqueAtomIndices,
                [len(np.where(result['equivalent_atoms'] == index)[0]) for index in uniqueAtomIndices]
                );

            # Store the symmetry tolerance used to call the spglib functions.

            self._symmetryAnalysisTolerance = tolerance;

    # -----------------------
    # Property Getter Methods
    # -----------------------

    def GetAtomCount(self):
        return len(self._atomData);

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

    def GetLatticeVectors(self):
        # Return the lattice vectors as a list of 1D NumPy arrays.

        return [np.copy(row) for row in self._latticeVectors[:]];

    def GetAtomTypeNumbers(self):
        # Return the atom-type numbers as a list of integers.

        return [atomTypeNumber for atomTypeNumber in self._atomTypeNumbersView[:]];

    def GetAtomPositions(self):
        # Return the atom positions as a list of 1D NumPy arrays.

        return [np.copy(row) for row in self._atomPositionsView[:]];

    def GetLatticeVectorsNumPy(self, copy = True):
        # Return the lattice vectors as a 3x3 NumPy array.

        if copy:
            return np.copy(self._latticeVectors);
        else:
            return self._latticeVectors;

    def GetAtomTypeNumbersNumPy(self, copy = True):
        # Return the atom-type numbers as a 1D NumPy array or the internal _atomTypeNumbersView field.

        if copy:
            return np.copy(self._atomTypeNumbersView);
        else:
            return self._atomTypeNumbersView;

    def GetAtomPositionsNumPy(self, copy = True):
        # Return the atom positions as an Nx3 NumPy array or the internal _atomPositionsView field.

        if copy:
            return np.copy(self._atomPositionsView);
        else:
            return self._atomPositionsView;

    def GetAtomDataNumPy(self, copy = True):
        # Return the internal _atomData field.

        if copy:
            return np.copy(self._atomData);
        else:
            return self._atomData;

    # -----------------------
    # Property Setter Methods
    # -----------------------

    def SetName(self, name):
        self._name = name;

    def SetAtom(self, index, atomType, atomPosition = None):
        # Sanity check.

        numAtoms = len(self._atomData);

        if index >= numAtoms:
            raise Exception("Error: Index {0} is out of range for number of atoms {1}.".format(index, numAtoms));

        # Convert to atomType number to an atom-type number.

        atomTypeNumber = AtomTypeToAtomTypeNumber(atomType);

        if atomTypeNumber != None:
            # Substitute atom.

            self._atomTypeNumbersView[index] = atomTypeNumber;

            if atomPosition != None:
                self._atomPositionsView[index] = [x % 1.0 for x in atomPosition];

            # Resort the atoms.

            self._atomData.sort();
        else:
            if atomPosition != None:
                warnings.warn("When atomType is set to None, the atom is removed and its position is not updated.", UserWarning);

            # Delete atom.

            numAtomsNew = numAtoms - 1;

            atomDataNew = np.zeros(numAtomsNew, dtype = Structure._AtomDataType);

            atomDataNew[:index] = self._atomData[:index];
            atomDataNew[index:] = self._atomData[index + 1:];

            # Update _atomData field and update views.

            self._atomData = atomDataNew;

            self._atomTypeNumbersView = atomDataNew.view(dtype = np.int64).reshape((numAtomsNew, 4))[:, 0];
            self._atomPositionsView = atomDataNew.view(dtype = np.float64).reshape((numAtomsNew, 4))[:, 1:];

        # Invalidate lazily-initialised property fields.

        self._ResetLazyPropertyFields();

    def SwapAtoms(self, atomType1, atomType2):
        atomTypeNumber1, atomTypeNumber2 = AtomTypeToAtomTypeNumber(atomType1), AtomTypeToAtomTypeNumber(atomType2);

        if atomType1 == None:
            raise Exception("Error: atomType1 cannot be None.");

        atomTypeNumbers = self._atomTypeNumbersView;

        if atomTypeNumber1 not in atomTypeNumbers:
            raise Exception("Error: atomType1 was not found in the structure.");

        if atomTypeNumber2 != None:
            # Using Boolean indexing prompts a copy, so we need to update the atom-type numbers in a loop.

            for index in np.argwhere(atomTypeNumbers == atomTypeNumber1):
                atomTypeNumbers[index] = atomTypeNumber2;

            # Updating the atom-type numbers may require a re-sort.

            self._atomData.sort();
        else:
            atomData = self._atomData;

            atomDataNew = atomData[atomTypeNumbers != atomType1];

            numAtomsNew = len(atomDataNew);

            # Since boolean indexing counts as "fancy indexing", we need to regenerate the type-number and positions views.

            self._atomData = atomDataNew;

            self._atomTypeNumbersView = atomDataNew.view(dtype = np.int64).reshape((numAtomsNew, 4))[:, 0];
            self._atomPositionsView = atomDataNew.view(dtype = np.float64).reshape((numAtomsNew, 4))[:, 1:];

    # ----------------
    # Symmetry Methods
    # ----------------

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

    # ------------------
    # Comparison Methods
    # ------------------

    def CompareLatticeVectors(self, structure, tolerance = None):
        # If no equivalence tolerance is supplied, use the default.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        return (np.abs(self._latticeVectors - structure.GetLatticeVectorsNumPy(copy = False)) <= tolerance).all();

    def CompareAtomTypeNumbers(self, structure):
        atomTypeNumbers1 = self._atomTypeNumbersView;
        atomTypeNumbers2 = structure.GetAtomTypeNumbersNumPy(copy = False);

        # First check that the number of atoms match.

        if len(atomTypeNumbers1) != len(atomTypeNumbers2):
            return False;

        # If they do, perform an element-wise comparison.

        return (atomTypeNumbers1 == atomTypeNumbers2).all();

    def CompareAtomPositions(self, structure, tolerance = None):
        # If no equivalence tolerance is supplied, use the default.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        atomPositions1 = self._atomPositionsView;
        atomPositions2 = structure.GetAtomPositionsNumPy(copy = False);

        # First check that the number of atoms match.

        if len(atomPositions1) != len(atomPositions2):
            return False;

        # If they do, perform an element-wise comparison using the tolerance.

        return (np.abs(atomPositions1 - atomPositions2) <= tolerance).all();

    def CompareStructure(self, structure, tolerance = None):
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

        # All tests passed -> the structures are equivalent.

        return True;

    # ---------------
    # Utility Methods
    # ---------------

    def Clone(self):
        # Given the way the constructor works, all the fields except the name (which is immutable anyway) should be deep copied.

        return Structure(
            self._latticeVectors, self._atomPositionsView, self._atomTypeNumbersView, name = self._name
            );

    def GetAtomPositionsCartesian(self):
        # Convert the atom positions to Cartesian coordinates and return as a list of NumPy arrays.

        if self._pAtomPositionsCartesian == None:
            # Load the lattice vectors.

            v1, v2, v3 = self._latticeVectors[:];

            # Convert the atom positions to Cartesian coordinates.

            atomPositionsCartesian = [];

            for x, y, z in self._atomPositionsView[:]:
                atomPositionsCartesian.append(
                    x * v1 + y * v2 + z * v3
                    );

            self._pAtomPositionsCartesian = atomPositionsCartesian;

        return self._pAtomPositionsCartesian;

    def GetAtomTypeNumberPlaceholder(self):
        # Return a negative atom-type number that isn't present among the atom-type numbers.

        return min(0, np.min(self._atomTypeNumbersView)) - 1;

    def GetAtomicSymbolsCounts(self, atomicSymbolLookupTable = None):
        # If atomicSymbolLookupTable is not set and the internal _pAtomSymbolsCounts has been initialised, we can simply return a copy of that.

        if atomicSymbolLookupTable == None and self._pAtomicSymbolsCounts != None:
            atomicSymbols, atomCounts = self._pAtomicSymbolsCounts;

            return (atomicSymbols[:], atomCounts[:]);

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

        # To avoid having to store a reference to the symbol lookup table, copies of the atom types/counts are only stored internally if atomicSymbolLookupTable is not set.

        if atomicSymbolLookupTable == None:
            self._pAtomicSymbolsCounts = (atomicSymbols[:], atomCounts[:]);

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

    def GetSupercell(self, supercellDim):
        dimA, dimB, dimC = supercellDim;

        # Sanity check.

        if dimA < 1 or dimB < 1 or dimC < 1:
            raise Exception("Error: Supercell dimensions must be >= 1.");

        if dimA % 1 != 0 or dimB % 1 != 0 or dimC % 1 != 0:
            raise Exception("Error: Supercell dimensions must be integers.");

        # If dimA = dimB = dimC = 1, we can simply return a clone.

        if dimA == 1 and dimB == 1 and dimC == 1:
            return self.Clone();

        atomTypeNumbers = self._atomTypeNumbersView;

        numAtoms = len(atomTypeNumbers);

        # Copy and rescale the atom positions.

        atomPositions = np.copy(self._atomPositionsView);

        atomPositions[:] /= np.array([dimA, dimB, dimC], dtype = np.float64);

        # Create arrays for the new atom type-numbers and positions.

        newAtomTypeNumbers = np.zeros(
            dimA * dimB * dimC * numAtoms, dtype = np.int64
            );

        newAtomPositions = np.zeros(
            (len(newAtomTypeNumbers), 3), dtype = np.float64
            );

        # Generate the new atom positions and list of atomic numbers for the supercell.

        pointer = 0;

        for i in range(0, dimA):
            for j in range(0, dimB):
                for k in range(0, dimC):
                    translationVector = np.array(
                        [float(i) / dimA, float(j) / dimB, float(k) / dimC], dtype = np.float64
                        );

                    newAtomTypeNumbers[pointer:pointer + numAtoms] = atomTypeNumbers;
                    newAtomPositions[pointer:pointer + numAtoms] = atomPositions + translationVector;

                    pointer += numAtoms;

        # Generate the new lattice vectors.

        newLatticeVectors = np.copy(self._latticeVectors);

        newLatticeVectors[0, :] *= dimA;
        newLatticeVectors[1, :] *= dimB;
        newLatticeVectors[2, :] *= dimC;

        # Generate a new name..

        newName = None;

        if self._name != None:
            # If the current structure has a name, generate a new one by appending "(<dimA>x<dimB>x<dimC> SC)" to it.

            newName = "{0} ({1}x{2}x{3} SC)".format(self._name, dimA, dimB, dimC);

        # Return a new structure.

        return Structure(
            newLatticeVectors, newAtomPositions, newAtomTypeNumbers, name = newName
            );

    # -------------
    # Static Fields
    # -------------

    _AtomDataType = [('TypeNumber', 'i8'), ('PosX', 'f8'), ('PosY', 'f8'), ('PosZ', 'f8')];

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
